"""Tests for cutadapt PCR-primer removal in the Preprocess stage."""

from __future__ import annotations

import gzip
import shutil
import tempfile
import subprocess
import unittest
from pathlib import Path

from core.models import AnalysisMode, Marker, Sample, SamplePaths
from preprocess import cutadapt_runner as cr

FWD = "CCAGCASCYGCGGTAATTCC"  # 18S V4 forward, with degenerate S and Y
REV = "ACTTTCGTTCTTGATYRA"  # 18S V4 reverse, with degenerate Y and R
# A real read carries concrete bases, never the ambiguity letters — those live only in
# the primer. These are one concrete instance of each primer (S->C, Y->C, R->A), so the
# reads are what an instrument would actually emit.
FWD_READ = "CCAGCACCCGCGGTAATTCC"
REV_READ = "ACTTTCGTTCTTGATCAA"
BODY = "ACGTACGTACGTACGTACGTACGTACGTAC"


class _Cfg:
    def __init__(self, out: Path, params: dict):
        self.output_dir = out
        self.vsearch_params = params


def _write_fastq(path: Path, seq: str, n: int = 5) -> None:
    with gzip.open(path, "wt") as handle:
        for i in range(n):
            handle.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def _sample(r1: Path, r2=None, mode=AnalysisMode.SE) -> Sample:
    return Sample("s1", "S1", SamplePaths(r1=r1, r2=r2), mode, Marker.OTHER)


class ResolvePrimersTests(unittest.TestCase):
    def test_explicit_fields_win_over_preset(self):
        cfg = _Cfg(Path("."), {"primer_fwd": "AAA", "primer_rev": "TTT", "primer_preset": "18S_V4_TAReuk"})
        self.assertEqual(cr.resolve_primers(cfg), ("AAA", "TTT"))

    def test_preset_used_when_no_explicit(self):
        cfg = _Cfg(Path("."), {"primer_preset": "18S_V4_TAReuk"})
        fwd, rev = cr.resolve_primers(cfg)
        self.assertTrue(fwd and rev)

    def test_nothing_configured(self):
        self.assertEqual(cr.resolve_primers(_Cfg(Path("."), {})), ("", ""))


class PassthroughTests(unittest.TestCase):
    def test_no_primers_returns_original_reads(self):
        with tempfile.TemporaryDirectory() as d:
            r1 = Path(d) / "r1.fastq.gz"
            _write_fastq(r1, FWD_READ + BODY)
            res = cr.run_cutadapt(_sample(r1), _Cfg(Path(d), {}))
            self.assertTrue(res.success)
            self.assertEqual(res.metrics["primers_removed"], 0.0)
            self.assertEqual(res.files["r1_primer"], r1)

    def test_cutadapt_absent_passes_through_with_warning(self):
        original = cr.cutadapt_available
        cr.cutadapt_available = lambda: False
        try:
            with tempfile.TemporaryDirectory() as d:
                r1 = Path(d) / "r1.fastq.gz"
                _write_fastq(r1, FWD_READ + BODY)
                res = cr.run_cutadapt(_sample(r1), _Cfg(Path(d), {"primer_fwd": FWD}))
            self.assertTrue(res.success)
            self.assertEqual(res.metrics["primers_removed"], 0.0)
            self.assertIn("cutadapt", res.error_message.lower())
        finally:
            cr.cutadapt_available = original


@unittest.skipIf(not cr.cutadapt_available(), "cutadapt not installed")
class RemovalTests(unittest.TestCase):
    def test_degenerate_primer_is_removed_from_every_read(self):
        with tempfile.TemporaryDirectory() as d:
            r1 = Path(d) / "r1.fastq.gz"
            _write_fastq(r1, FWD_READ + BODY, n=5)
            res = cr.run_cutadapt(_sample(r1), _Cfg(Path(d), {"primer_fwd": FWD, "threads": 2}))
            self.assertEqual(res.metrics["primers_removed"], 1.0)
            self.assertEqual(res.metrics["reads_with_primer"], 5.0)
            out = Path(res.files["r1_primer"])
            trimmed = gzip.open(out, "rt").readlines()[1].strip()
            self.assertEqual(trimmed, BODY)

    def test_pe_removes_from_both_reads(self):
        with tempfile.TemporaryDirectory() as d:
            r1 = Path(d) / "r1.fastq.gz"
            r2 = Path(d) / "r2.fastq.gz"
            _write_fastq(r1, FWD_READ + BODY)
            _write_fastq(r2, REV_READ + BODY)
            res = cr.run_cutadapt(_sample(r1, r2, AnalysisMode.PE), _Cfg(Path(d), {"primer_fwd": FWD, "primer_rev": REV}))
            self.assertEqual(res.metrics["primers_removed"], 1.0)
            self.assertTrue(Path(res.files["r1_primer"]).exists())
            self.assertTrue(Path(res.files["r2_primer"]).exists())
            self.assertEqual(gzip.open(res.files["r1_primer"], "rt").readlines()[1].strip(), BODY)
            self.assertEqual(gzip.open(res.files["r2_primer"], "rt").readlines()[1].strip(), BODY)


class NoDoubleStripTests(unittest.TestCase):
    def test_vsearch_skips_length_strip_when_primers_removed_upstream(self):
        from analysis.vsearch_runner import _primer_strip_lengths

        params = {"primer_fwd": FWD, "primer_rev": "ACTTTCGTTCTTGATYRA"}
        self.assertNotEqual(_primer_strip_lengths(params), (0, 0))
        params["primers_removed_upstream"] = True
        self.assertEqual(_primer_strip_lengths(params), (0, 0))


if __name__ == "__main__":
    unittest.main()


class PrimerPresetDefaultTests(unittest.TestCase):
    """The default preset was "" — the same value as "do not strip primers".

    So a run that never touched the control was indistinguishable from one that opted
    out, and got no primer removal: 7/7 OTUs on run 20260716T080902Z still carried
    TAReuk454FWD1 plus the heterogeneous 5' spacer in front of it, which split one cod
    (Micromesistius poutassou, >99% on both) across OTU1 and OTU2.
    """

    def _config(self, marker, **params):
        from core.models import AnalysisBranch, AnalysisMode, Marker, RunConfig, Sample, SamplePaths

        sample = Sample(
            internal_id="s1",
            display_label="s1",
            paths=SamplePaths(r1=Path("/tmp/r1.fastq")),
            mode=AnalysisMode.SE,
            marker=marker,
        )
        config = RunConfig(run_id="r", samples=[sample], branch=AnalysisBranch.VSEARCH)
        config.vsearch_params.update(params)
        return config

    def test_auto_resolves_the_tareuk_pair_for_an_18s_run(self):
        from core.models import Marker
        from preprocess.cutadapt_runner import resolve_primers

        forward, reverse = resolve_primers(self._config(Marker.SSU_18S, primer_preset="auto"))
        self.assertEqual(forward, "CCAGCASCYGCGGTAATTCC")
        self.assertEqual(reverse, "ACTTTCGTTCTTGATYRA")

    def test_an_absent_preset_is_treated_as_auto(self):
        # Runs stored before "auto" existed have no primer_preset key at all.
        from core.models import Marker
        from preprocess.cutadapt_runner import resolve_primers

        config = self._config(Marker.SSU_18S)
        config.vsearch_params.pop("primer_preset", None)
        self.assertEqual(resolve_primers(config)[0], "CCAGCASCYGCGGTAATTCC")

    def test_an_explicit_none_is_obeyed_and_not_auto_filled(self):
        # "" is the user choosing not to strip primers. Overriding that would be worse
        # than the bug this fixes.
        from core.models import Marker
        from preprocess.cutadapt_runner import resolve_primers

        self.assertEqual(resolve_primers(self._config(Marker.SSU_18S, primer_preset="")), ("", ""))

    def test_coi_has_no_auto_default(self):
        # Folmer and Leray are both standard and amplify different regions; there is no
        # single pair to fall back on, so auto declines rather than picks one.
        from core.models import Marker
        from preprocess.cutadapt_runner import resolve_primers

        self.assertEqual(resolve_primers(self._config(Marker.COI, primer_preset="auto")), ("", ""))

    def test_explicit_sequences_win_over_auto(self):
        from core.models import Marker
        from preprocess.cutadapt_runner import resolve_primers

        config = self._config(Marker.SSU_18S, primer_preset="auto", primer_fwd="ACGT", primer_rev="TTTT")
        self.assertEqual(resolve_primers(config), ("ACGT", "TTTT"))

    def test_a_mixed_marker_run_gets_no_auto_pair(self):
        from core.models import AnalysisBranch, AnalysisMode, Marker, RunConfig, Sample, SamplePaths
        from preprocess.cutadapt_runner import resolve_primers

        samples = [
            Sample(internal_id=f"s{i}", display_label=f"s{i}", paths=SamplePaths(r1=Path("/tmp/r.fastq")), mode=AnalysisMode.SE, marker=m)
            for i, m in enumerate((Marker.SSU_18S, Marker.COI))
        ]
        config = RunConfig(run_id="r", samples=samples, branch=AnalysisBranch.VSEARCH)
        config.vsearch_params["primer_preset"] = "auto"
        self.assertEqual(resolve_primers(config), ("", ""))

    def test_the_setup_default_is_auto_not_the_do_not_strip_value(self):
        """Guards the actual regression: the dropdown shipped with value=""."""
        import re

        source = Path("ui/pages/page_setup.py").read_text(encoding="utf-8")
        block = source[source.index('id="setup-primer-preset"'):]
        value = re.search(r'value="([^"]*)"', block).group(1)
        self.assertEqual(value, "auto")


class DiscardUntrimmedTests(unittest.TestCase):
    """Reads with no primer match are dropped, not passed through uncut.

    Provenance, not signal: a read cutadapt could not trim cannot be shown to be an
    amplicon of the targeted region. Measured on the COI/Leray run, 92.32% of input
    carried the forward primer, so this discards ~7.7% of INPUT — but only 0.19% of the
    reads that currently reach clustering, because assembly and quality filtering
    already remove most of the rest.
    """

    def _command(self, mode, r2=None):
        import tempfile
        from unittest.mock import patch

        from core.models import AnalysisMode, AnalysisBranch, RunConfig, Sample, SamplePaths
        from preprocess.cutadapt_runner import run_cutadapt

        with tempfile.TemporaryDirectory() as d:
            from pathlib import Path

            r1 = Path(d) / "s_R1.fastq"
            r1.write_text("@a\nACGT\n+\nIIII\n", encoding="utf-8")
            r2_path = None
            if r2:
                r2_path = Path(d) / "s_R2.fastq"
                r2_path.write_text("@a\nACGT\n+\nIIII\n", encoding="utf-8")
            sample = Sample(
                internal_id="s", display_label="S",
                paths=SamplePaths(r1=r1, r2=r2_path), mode=mode, marker=Marker.COI,
            )
            config = RunConfig(run_id="t", samples=[sample], branch=AnalysisBranch.VSEARCH)
            config.output_dir = Path(d) / "out"
            config.vsearch_params["primer_preset"] = "COI_Leray"
            captured = {}

            def fake_run(command, **kwargs):
                captured["command"] = command
                return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

            with patch("preprocess.cutadapt_runner.cutadapt_available", return_value=True), \
                patch("preprocess.cutadapt_runner.subprocess.run", side_effect=fake_run):
                run_cutadapt(sample, config)
            return captured.get("command", [])

    def test_paired_end_discards_untrimmed_pairs(self):
        from core.models import AnalysisMode

        command = self._command(AnalysisMode.PE, r2=True)

        self.assertIn("--discard-untrimmed", command)
        self.assertIn("--pair-filter=any", command)

    def test_single_end_discards_untrimmed_reads_without_pair_filter(self):
        """cutadapt REJECTS --pair-filter in single-end mode. Passing it there made the
        step fail and fall back to untrimmed passthrough — primers silently kept."""
        from core.models import AnalysisMode

        command = self._command(AnalysisMode.SE)
        self.assertIn("--discard-untrimmed", command)
        self.assertNotIn("--pair-filter=any", command)

    def test_the_flag_is_a_configurable_choice_not_a_literal(self):
        """A scientific choice is a parameter. Same lesson as MIN_UNIQUE_SIZE — which
        was decided as 1 and then silently overridden by a 2 typed into the UI."""
        from unittest.mock import patch

        from core.models import AnalysisMode

        with patch("preprocess.cutadapt_runner.DISCARD_UNTRIMMED", False):
            self.assertNotIn("--discard-untrimmed", self._command(AnalysisMode.SE))

    def test_default_is_on(self):
        from core.config import DISCARD_UNTRIMMED

        self.assertTrue(DISCARD_UNTRIMMED)
