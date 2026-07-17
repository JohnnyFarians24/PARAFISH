"""Marker-filtered reference tests.

These cover the ways a reference filter can be quietly *wrong* rather than broken: a
reference that still classifies, but against the wrong region, the wrong strand or a
corrupted taxonomy would produce confident nonsense instead of an error.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from core.preflight import preflight_notes
from databases.reference_filter import (
    FILTER_FORMAT_VERSION,
    build_marker_reference,
    marker_filter_applies,
    normalize_marker,
    resolve_marker,
    reverse_complement,
)

FORWARD = "CCAGCASCYGCGGTAATTCC"
REVERSE = "ACTTTCGTTCTTGATYRA"

# A synthetic SSU-like record: filler, the V4 flanks, and the region between them.
V4_INSERT = "ACGTTGCATGCCTAGGTACCAGATTCCGGATTACGCATGCATCGGATCCAAGTTACGCA" * 4
LEADING = "TTACGGTCAGCATTAGCCTAGGACT" * 3
TRAILING = "GGATCCATTGCACTGGATCAGGTTA" * 3


def _record(seq_id: str, lineage: str, sequence: str) -> str:
    return f">{seq_id};tax={lineage};\n{sequence}\n"


class ReverseComplementTests(unittest.TestCase):
    def test_iupac_degenerate_bases_are_complemented(self):
        # The reverse primer is ordered 5'->3' on the opposite strand; getting R/Y or
        # K/M backwards yields a primer that silently matches nothing.
        self.assertEqual(reverse_complement(REVERSE), "TYRATCAAGAACGAAAGT")
        self.assertEqual(reverse_complement("RYSWKMBDHVN"), "NBDHVKMWSRY")

    def test_reverse_complement_is_its_own_inverse(self):
        for seq in (FORWARD, REVERSE, "ACGT", "RYSWKMBDHVN"):
            self.assertEqual(reverse_complement(reverse_complement(seq)), seq)

    def test_rna_uracil_complements_to_adenine_and_comes_back_as_dna(self):
        # SILVA ships RNA; U must behave as T or the primer match fails outright. The
        # round trip is deliberately asymmetric (U -> A -> T): everything downstream of
        # the filter is DNA, so there is nothing to convert back to.
        self.assertEqual(reverse_complement("UUUU"), "AAAA")
        self.assertEqual(reverse_complement(reverse_complement("ACGTU")), "ACGTT")


class MarkerResolutionTests(unittest.TestCase):
    def test_marker_labels_normalize_to_scope_keys(self):
        self.assertEqual(normalize_marker("18S"), "18S")
        self.assertEqual(normalize_marker("SSU_18S"), "18S")
        self.assertEqual(normalize_marker(""), "")

    def test_mixed_marker_run_resolves_to_nothing(self):
        # No single amplicon covers both, so such a run must get the unfiltered
        # reference rather than one that is wrong for half its samples.
        self.assertEqual(resolve_marker(["18S", "COI"]), "")
        self.assertEqual(resolve_marker(["18S", "18S"]), "18S")
        self.assertEqual(resolve_marker([]), "")

    def test_coi_without_primers_gets_no_filter(self):
        # SILVA/PR2 hold no COI at all, and there is no lineage scope to apply.
        self.assertFalse(marker_filter_applies("COI", ("", "")))
        self.assertTrue(marker_filter_applies("18S", ("", "")))


class BuildMarkerReferenceTests(unittest.TestCase):
    def _silva_like(self, path: Path) -> Path:
        """A mixed-domain reference in RNA, as SILVA ships it."""
        records = [
            # Taxon name deliberately contains a U: rewriting U->T in headers would
            # corrupt this into 'Tnruhdinium'.
            _record("EUK1", "d:Eukaryota,p:Alveolata,g:Unruhdinium,s:Unruhdinium_kevei", (LEADING + FORWARD + V4_INSERT + reverse_complement(REVERSE) + TRAILING).replace("T", "U")),
            _record("EUK2", "d:Eukaryota,p:Amorphea,g:Teleostei,s:Gadus_morhua", (LEADING + FORWARD + V4_INSERT + reverse_complement(REVERSE) + TRAILING).replace("T", "U")),
            _record("BAC1", "d:Bacteria,p:Proteobacteria,s:Escherichia_coli", (LEADING + FORWARD + V4_INSERT + reverse_complement(REVERSE) + TRAILING).replace("T", "U")),
            _record("ARC1", "d:Archaea,p:Euryarchaeota,s:Methanobrevibacter_smithii", (LEADING + V4_INSERT + TRAILING).replace("T", "U")),
        ]
        path.write_text("".join(records), encoding="utf-8")
        return path

    def test_filter_keeps_eukaryotes_cuts_amplicon_and_preserves_taxonomy(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = self._silva_like(Path(temp_dir) / "ref.sintax.fasta")
            reference = build_marker_reference(source, "18S", (FORWARD, REVERSE))

            self.assertIsNotNone(reference)
            text = reference.path.read_text(encoding="utf-8")
            self.assertIn("EUK1", text)
            self.assertIn("EUK2", text)
            # Bacteria/archaea are 89% of SILVA and 0% of an 18S answer.
            self.assertNotIn("BAC1", text)
            self.assertNotIn("ARC1", text)
            # The U in the taxon name survived; the RNA in the sequence did not.
            self.assertIn("Unruhdinium_kevei", text)
            sequences = [line for line in text.splitlines() if not line.startswith(">")]
            self.assertTrue(all("U" not in seq for seq in sequences))
            # The primers are cut off with the flanks: what is left is the region
            # between them, which is the region the reads actually span.
            self.assertTrue(all(seq == V4_INSERT for seq in sequences))
            self.assertLess(reference.size_bytes, source.stat().st_size)
            self.assertEqual(reference.records, 2)
            self.assertEqual(reference.source_records, 4)

    def test_reference_is_cached_and_carries_the_format_version(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = self._silva_like(Path(temp_dir) / "ref.sintax.fasta")
            first = build_marker_reference(source, "18S", (FORWARD, REVERSE))
            self.assertIn(f".filtered.{FILTER_FORMAT_VERSION}.fasta", first.path.name)

            # A cache hit must not shell out to cutadapt again.
            with patch("databases.reference_filter._extract_amplicon") as extract:
                second = build_marker_reference(source, "18S", (FORWARD, REVERSE))
            extract.assert_not_called()
            self.assertEqual(second.path, first.path)
            self.assertEqual(second.records, first.records)

    def test_a_different_primer_pair_gets_a_different_cache(self):
        # A reference cut for V4 is not the reference a V9 run needs; sharing one cache
        # entry between them would hand a run the wrong region without a word.
        with tempfile.TemporaryDirectory() as temp_dir:
            source = self._silva_like(Path(temp_dir) / "ref.sintax.fasta")
            v4 = build_marker_reference(source, "18S", (FORWARD, REVERSE))
            v9 = build_marker_reference(source, "18S", ("GTACACACCGCCCGTC", "TGATCCTTCTGCAGGTTCACCTAC"))
            self.assertNotEqual(v4.path, v9.path)

    def test_stale_cache_is_rebuilt_when_the_source_changes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = self._silva_like(Path(temp_dir) / "ref.sintax.fasta")
            first = build_marker_reference(source, "18S", (FORWARD, REVERSE))
            self.assertEqual(first.records, 2)

            extra = _record("EUK3", "d:Eukaryota,p:Rhizaria,s:Bigelowiella_natans", LEADING + FORWARD + V4_INSERT + reverse_complement(REVERSE) + TRAILING)
            source.write_text(source.read_text(encoding="utf-8") + extra, encoding="utf-8")
            import os

            os.utime(source, (first.path.stat().st_mtime + 10, first.path.stat().st_mtime + 10))
            rebuilt = build_marker_reference(source, "18S", (FORWARD, REVERSE))
            self.assertEqual(rebuilt.records, 3)

    def test_primers_that_match_nothing_fall_back_to_the_lineage_filter(self):
        # A primer pair from another marker would otherwise produce a correctly-sized
        # reference that classifies nothing at all.
        with tempfile.TemporaryDirectory() as temp_dir:
            source = self._silva_like(Path(temp_dir) / "ref.sintax.fasta")
            reference = build_marker_reference(source, "18S", ("GGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCC"))
            self.assertIsNotNone(reference)
            text = reference.path.read_text(encoding="utf-8")
            self.assertIn("EUK1", text)
            self.assertNotIn("BAC1", text)
            # Rung 2 of the ladder: still eukaryote-only, but sequences kept whole.
            self.assertEqual(reference.records, 2)
            self.assertTrue(any("primers" in step for step in reference.steps))

    def test_all_eukaryote_reference_without_primers_is_left_alone(self):
        # PR2 is already eukaryote-only: "filtering" it without primers would write a
        # full-size copy to disk and call it filtered.
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "pr2.sintax.fasta"
            source.write_text(
                _record("P1", "d:Eukaryota,p:Alveolata,s:Unruhdinium_kevei", V4_INSERT)
                + _record("P2", "d:Eukaryota,p:Amorphea,s:Gadus_morhua", V4_INSERT),
                encoding="utf-8",
            )
            self.assertIsNone(build_marker_reference(source, "18S", ("", "")))

    def test_missing_cutadapt_still_yields_the_lineage_filter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = self._silva_like(Path(temp_dir) / "ref.sintax.fasta")
            with patch("databases.reference_filter.cutadapt_available", return_value=False):
                reference = build_marker_reference(source, "18S", (FORWARD, REVERSE))
            self.assertIsNotNone(reference)
            self.assertEqual(reference.records, 2)
            self.assertNotIn("BAC1", reference.path.read_text(encoding="utf-8"))

    def test_unknown_marker_gets_no_filter(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = self._silva_like(Path(temp_dir) / "ref.sintax.fasta")
            self.assertIsNone(build_marker_reference(source, "", ("", "")))
            self.assertIsNone(build_marker_reference(Path(temp_dir) / "absent.fasta", "18S", (FORWARD, REVERSE)))


class PreflightMarkerNoteTests(unittest.TestCase):
    class _Paths:
        pr2 = None
        blast = None
        bold_dir = None

        def __init__(self, silva: Path):
            self.silva = silva

    def _notes(self, temp_dir: str, marker_filter: str) -> list[str]:
        source = Path(temp_dir) / "silva.fasta"
        source.write_text("A" * 4096, encoding="utf-8")
        with patch("core.preflight.available_memory_bytes", return_value=1024):
            return preflight_notes(["silva"], self._Paths(source), 1, marker_filter=marker_filter)

    def test_oversized_reference_warns_when_no_filter_applies(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            notes = " ".join(self._notes(temp_dir, ""))
            self.assertIn("ignorada", notes)

    def test_no_skip_warning_once_the_reference_will_be_filtered(self):
        # The full size is real but no longer the size that gets loaded; warning about
        # it would send the user to fix a problem that no longer exists.
        with tempfile.TemporaryDirectory() as temp_dir:
            notes = " ".join(self._notes(temp_dir, "18S"))
            self.assertNotIn("ignorada", notes)
            self.assertIn("18S", notes)
            self.assertIn("reduzida", notes)


if __name__ == "__main__":
    unittest.main()
