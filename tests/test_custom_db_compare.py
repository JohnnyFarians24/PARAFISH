"""Tests for the custom-database BLAST comparison feature."""

from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

from analysis.custom_db_compare import NO_MATCH_LABEL, build_custom_comparison, write_detection_table
from databases.custom_db import MIN_REFERENCE_LENGTH, CustomDBConfig, blast_tools_available, detect_fasta_kind, reference_label, run_custom_blast


class DetectKindTests(unittest.TestCase):
    # These fixtures used a 4 bp "ACGT" as a stand-in body while detect_fasta_kind read
    # only headers. Sequence validation now applies, and a 4 bp subject is unusable by
    # blastn (word size >= 11) — so the placeholder is grown to something a user could
    # actually upload. The assertions below are unchanged.
    BODY = "ACGT" * 20

    def _write(self, path: Path, text: str) -> Path:
        path.write_text(text, encoding="utf-8")
        return path

    def test_plain_fasta_is_blast_kind(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write(Path(d) / "ref.fasta", f">seq1 Genus species\n{self.BODY}\n>seq2 Other taxon\n{self.BODY}\n")
            kind, error = detect_fasta_kind(p)
            self.assertEqual(kind, "blast")
            self.assertEqual(error, "")

    def test_sintax_fasta_is_sintax_kind(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write(Path(d) / "ref.fasta", f">seq1;tax=d:Eukaryota,g:Diplostomum,s:Diplostomum spathaceum;\n{self.BODY}\n")
            kind, error = detect_fasta_kind(p)
            self.assertEqual(kind, "sintax")

    def test_non_fasta_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write(Path(d) / "ref.txt", "not a fasta\n")
            kind, error = detect_fasta_kind(p)
            self.assertEqual(kind, "")
            self.assertTrue(error)

    def test_config_validate_accepts_plain_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write(Path(d) / "ref.fasta", f">seq1 Genus species\n{self.BODY}\n")
            cfg = CustomDBConfig(p)
            valid, error = cfg.validate()
            self.assertTrue(valid, error)
            self.assertEqual(cfg.kind, "blast")


class ReferenceLabelTests(unittest.TestCase):
    def test_sintax_species_field(self):
        self.assertEqual(reference_label("ID1;tax=d:Eukaryota,g:Diplostomum,s:Diplostomum_spathaceum;"), "Diplostomum spathaceum")

    def test_ncbi_style_genus_species(self):
        self.assertEqual(reference_label("ACC123 Gyrodactylus salaris voucher X"), "Gyrodactylus salaris")

    def test_bare_id(self):
        self.assertEqual(reference_label("ACC123"), "ACC123")


class BuildComparisonTests(unittest.TestCase):
    def _setup(self, d: Path):
        besthit = d / "besthit.csv"
        besthit.write_text(
            "OTU_ID,subject_id,species,pident,length,evalue,bitscore,status\n"
            "OTU1,REF1,Diplostomum spathaceum,99.5,300,0.0,550,Match\n"
            "OTU2,REF2,Gyrodactylus salaris,88.0,290,1e-40,300,Match\n"
            "OTU3,,None,,,,,None\n",
            encoding="utf-8",
        )
        matrix = d / "otu_matrix.csv"
        matrix.write_text(
            "#OTU ID,S1,S2\n"
            "OTU1,100,0\n"
            "OTU2,10,40\n"
            "OTU3,0,25\n",
            encoding="utf-8",
        )
        return besthit, matrix

    def test_detection_matrix_and_summary(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            comp = build_custom_comparison(besthit, matrix)
            self.assertNotIn("error", comp)
            self.assertEqual(comp["summary"]["n_otus"], 3)
            self.assertEqual(comp["summary"]["n_matched"], 2)
            self.assertEqual(comp["summary"]["n_unmatched"], 1)
            self.assertEqual(comp["samples"], ["S1", "S2"])
            self.assertIn(NO_MATCH_LABEL, comp["taxa"])
            # No match row holds OTU3's reads (25 in S2).
            det = {name: row for name, row in zip(comp["detection"]["taxa"], comp["detection"]["values"])}
            self.assertEqual(det[NO_MATCH_LABEL], [0, 25])
            self.assertEqual(det["Diplostomum spathaceum"], [100, 0])

    def test_min_identity_demotes_weak_hits(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            comp = build_custom_comparison(besthit, matrix, min_identity=95.0)
            # OTU2 at 88% now counts as unmatched.
            self.assertEqual(comp["summary"]["n_matched"], 1)
            otu2 = next(r for r in comp["per_otu"] if r["OTU_ID"] == "OTU2")
            self.assertEqual(otu2["status"], "None")

    def test_species_filter(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            comp = build_custom_comparison(besthit, matrix, species_filter="diplostomum")
            self.assertEqual([r["OTU_ID"] for r in comp["per_otu"]], ["OTU1"])
            self.assertEqual(comp["taxa"], ["Diplostomum spathaceum"])

    def test_public_species_from_consensus(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            consensus = d / "consensus.csv"
            consensus.write_text(
                "ASV_ID,Consensus_Species\nOTU1,Diplostomum spathaceum\nOTU2,\n",
                encoding="utf-8",
            )
            comp = build_custom_comparison(besthit, matrix, consensus_csv=consensus)
            by_id = {r["OTU_ID"]: r for r in comp["per_otu"]}
            self.assertEqual(by_id["OTU1"]["public_species"], "Diplostomum spathaceum")
            self.assertEqual(by_id["OTU2"]["public_species"], "None")

    def test_write_detection_table(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            comp = build_custom_comparison(besthit, matrix)
            out = write_detection_table(comp["detection"], d / "detection.csv")
            self.assertTrue(out.exists())
            text = out.read_text(encoding="utf-8")
            self.assertIn("taxon,S1,S2", text)
            self.assertIn("No match", text)

    def test_missing_besthit_returns_error(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            _, matrix = self._setup(d)
            comp = build_custom_comparison(d / "nope.csv", matrix)
            self.assertIn("error", comp)


@unittest.skipUnless(blast_tools_available(), "blastn/makeblastdb not available")
class RealBlastTests(unittest.TestCase):
    def test_blast_reports_none_for_unmatched_query(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            ref_seq = (
                "TCGACCTGGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCAA"
                "TTATACCGCGAATGGCTCATTATAACAGTCGCACGTACCTTGACTATCTGGATAACCGTG"
                "GTAATTCTAGAGCTAATACATGCCGACGAGCGCTGACCCCCGGGGATGCGTGCATTTATC"
            )
            other_seq = (
                "GAACCTTACCAGGTCCAGACATAGTAAGGATTGACAGATTGAGAGCTCTTTCTTGATTTT"
                "GTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGAT"
                "AACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCATTCGCTAGTTGGCATCG"
            )
            reference = d / "ref.fasta"
            reference.write_text(f">REF1 Testus specius\n{ref_seq}\n", encoding="utf-8")
            query = d / "query.fasta"
            # OTU1 == reference (matches), OTU2 is unrelated (no hit).
            query.write_text(
                f">OTU1;size=5\n{ref_seq}\n>OTU2;size=3\n{other_seq}\n",
                encoding="utf-8",
            )
            result = run_custom_blast(query, reference, d / "out", perc_identity=90.0, threads=1)
            self.assertTrue(result.success, result.error_message)
            self.assertEqual(result.metrics["matched"], 1.0)
            self.assertEqual(result.metrics["total"], 2.0)
            besthit = Path(result.files["custom_blast_besthit_csv"])
            text = besthit.read_text(encoding="utf-8")
            self.assertIn("OTU1", text)
            self.assertIn("Testus specius", text)
            # The unmatched OTU2 row is present with a None status.
            import csv as _csv

            rows = {r["OTU_ID"]: r for r in _csv.DictReader(text.splitlines())}
            self.assertEqual(rows["OTU2"]["status"], "None")
            self.assertEqual(rows["OTU2"]["species"], "None")


if __name__ == "__main__":
    unittest.main()


class SequenceValidationTests(unittest.TestCase):
    """Task 3(a): the upload validates SEQUENCES, not just headers.

    detect_fasta_kind only ever read lines starting with ">", so a file with valid
    headers and protein, empty or junk bodies passed the upload and blew up a whole
    stage later inside makeblastdb — where the user can no longer fix it.
    """

    def _kind(self, text: str):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "ref.fasta"
            path.write_text(text, encoding="utf-8")
            return detect_fasta_kind(path)

    OK = "ACGT" * 20

    def test_protein_fasta_is_rejected_at_upload(self):
        kind, error = self._kind(">p1 some protein\nMKWVTFISLLLLFSSAYSRGV\n")

        self.assertEqual(kind, "")
        self.assertIn("non-nucleotide", error)

    def test_empty_sequence_body_is_rejected(self):
        kind, error = self._kind(f">s1\n\n>s2\n{self.OK}\n")

        self.assertEqual(kind, "")
        self.assertIn("has no bases", error)

    def test_headers_with_no_sequence_at_all_are_rejected(self):
        kind, error = self._kind(">s1\n>s2\n")

        self.assertEqual(kind, "")
        self.assertIn("has no bases", error)

    def test_junk_characters_are_rejected(self):
        kind, error = self._kind(">s1\n!!!???###\n")

        self.assertEqual(kind, "")
        self.assertIn("non-nucleotide", error)

    def test_a_sequence_below_the_floor_is_rejected(self):
        kind, error = self._kind(">s1\nACGTACGT\n")

        self.assertEqual(kind, "")
        self.assertIn("minimum", error)

    def test_the_length_floor_value_itself_is_pinned(self):
        """Written separately from the boundary test below, and deliberately literal.
        The boundary test derives its inputs FROM the constant, so it proves the
        comparison is < and not <=, but would pass just as happily with the floor at
        100. This one pins the number: 50 is the smallest reference worth keeping when
        blastn seeds on 11-base words and megablast on 28. Changing it is a scientific
        decision and should have to change a test."""
        self.assertEqual(MIN_REFERENCE_LENGTH, 50)

    def test_the_length_floor_is_pinned_at_its_exact_boundary(self):
        """Without this, MIN_REFERENCE_LENGTH could be 10 or 100 and the test above
        would still pass — the floor would have no proven value and the < / <= choice
        would be unproven. 50 is the smallest reference blastn can seed on usefully."""
        just_under = ("ACGTACGTAC" * 10)[: MIN_REFERENCE_LENGTH - 1]
        exactly = ("ACGTACGTAC" * 10)[:MIN_REFERENCE_LENGTH]

        self.assertEqual(self._kind(f">s1\n{just_under}\n")[0], "")
        self.assertEqual(self._kind(f">s1\n{exactly}\n"), ("blast", ""))

    def test_gaps_do_not_count_towards_the_length_floor(self):
        """Gaps are alignment padding, not bases. Padding a 49-base record with dashes
        must not buy it past a floor that exists to guarantee alignable sequence."""
        short = ("ACGTACGTAC" * 10)[: MIN_REFERENCE_LENGTH - 1]

        self.assertEqual(self._kind(f">s1\n{short}-----\n")[0], "")
        self.assertEqual(self._kind(f">s1\n{short}A-----\n")[0], "blast")

    def test_a_peptide_of_only_iupac_letters_is_caught_by_composition(self):
        """K, W, V, S, Y, R, M are all valid IUPAC codes — the alphabet check alone
        would wave this through."""
        kind, error = self._kind(">s1\n" + "RYKMSWBDHV" * 8 + "\n")

        self.assertEqual(kind, "")
        self.assertIn("A/C/G/T", error)

    def test_valid_dna_still_passes(self):
        self.assertEqual(self._kind(f">s1\n{self.OK}\n"), ("blast", ""))

    def test_wrapped_sequences_are_joined_before_the_length_check(self):
        """60-column FASTA is the norm; measuring line by line would reject all of it."""
        wrapped = ">s1\n" + "\n".join("ACGT" * 5 for _ in range(4)) + "\n"
        self.assertEqual(self._kind(wrapped), ("blast", ""))

    def test_ambiguity_codes_and_gaps_are_accepted(self):
        """A masked or aligned reference is still DNA. The validator must not overreach."""
        self.assertEqual(self._kind(">s1\n" + "ACGTN" * 15 + "---\n")[0], "blast")

    def test_sintax_reference_is_validated_too(self):
        kind, _ = self._kind(f">s1;tax=d:Eukaryota,g:Gadus,s:Gadus_morhua;\n{self.OK}\n")
        self.assertEqual(kind, "sintax")
        self.assertEqual(self._kind(">s1;tax=d:Eukaryota;\nMKWVTFISLLLL\n")[0], "")


class ZeroByteIsNotFailureTests(unittest.TestCase):
    """Task 3(b): with BLAST, zero hits legitimately writes a zero-byte table.

    Observed on a real run: custom_blast_besthit.csv was complete and correct, every
    feature honestly marked "None", and source_status said `failed` — because the
    reconciliation judged by the raw table's SIZE. The evidence for failure lives in
    the log, where blastn prints "(Severe Error) ... search failed" while exiting 0.
    """

    def _status(self, raw_text: str, log_text: str):
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "custom_blast_raw.tsv"
            raw.write_text(raw_text, encoding="utf-8")
            log = Path(d) / "custom_blast.log"
            log.write_text(log_text, encoding="utf-8")
            entry = {"source": "Base customizada (BLAST)", "status": "no_hits",
                     "detail": "0/112", "key": "CUSTOM", "raw": raw, "log": log}
            return _reconcile_source_status([entry], {}, 112, marker="COI")[0]

    def test_empty_output_with_a_clean_log_is_no_hits(self):
        out = self._status("", "Custom BLAST: 0/112 OTUs com match.\n")

        self.assertEqual(out["status"], "no_hits")

    def test_empty_output_with_a_blast_error_in_the_log_is_failed(self):
        out = self._status("", "Error: (Severe Error) Blast search error: search failed\n")

        self.assertEqual(out["status"], "failed")
        self.assertIn("search failed", out["detail"])

    def test_a_populated_table_is_still_no_hits_when_nothing_was_assigned(self):
        row = "OTU1\tsub\t99\t300\t0\t0\t1\t300\t1\t300\t1e-50\t500\t100\n"
        self.assertEqual(self._status(row, "ok\n")["status"], "no_hits")

    def test_sources_without_a_log_keep_the_size_heuristic(self):
        """NCBI and BOLD pass no log here; their existing behaviour must not move."""
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "bold.parquet"
            raw.write_text("", encoding="utf-8")
            entry = {"source": "BOLD", "status": "ok", "detail": "x", "key": "BOLD", "raw": raw}
            self.assertEqual(_reconcile_source_status([entry], {}, 7, marker="COI")[0]["status"], "failed")


class ExtensionGatesAgreeTests(unittest.TestCase):
    """Task 3(c): three gates guarded one upload and disagreed about .fna."""

    def test_fna_is_accepted_by_every_gate(self):
        from ui.callbacks_databases import _is_supported_custom_db_filename

        for name in ("ref.fna", "ref.fna.gz", "ref.fasta", "ref.fa", "ref.fasta.gz", "ref.fa.gz"):
            with self.subTest(name=name):
                self.assertTrue(_is_supported_custom_db_filename(name))

    def test_the_validator_agrees_with_the_ui_gate(self):
        from ui.callbacks_databases import SUPPORTED_CUSTOM_DB_SUFFIXES

        with tempfile.TemporaryDirectory() as d:
            for suffix in SUPPORTED_CUSTOM_DB_SUFFIXES:
                if suffix.endswith(".gz"):
                    continue
                path = Path(d) / f"ref{suffix}"
                path.write_text(">s1\n" + "ACGT" * 20 + "\n", encoding="utf-8")
                ok, error = CustomDBConfig(path).validate()
                with self.subTest(suffix=suffix):
                    self.assertTrue(ok, error)

    def test_an_unsupported_extension_is_still_refused(self):
        from ui.callbacks_databases import _is_supported_custom_db_filename

        self.assertFalse(_is_supported_custom_db_filename("ref.txt"))
        self.assertFalse(_is_supported_custom_db_filename("ref.xlsx"))
