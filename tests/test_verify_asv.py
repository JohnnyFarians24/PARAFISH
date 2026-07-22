"""verify.py must SEE an ASV run — the instrument that verifies the instrument, seen to bite.

verify.py was OTU-blind: every check read ``databases/otu/`` and ``analysis/pooled/otu_table.tsv``,
so on an ASV run each SKIPped and the run exited 0 — a blind instrument goes quiet, it does not
scream, so an ASV run's "verify passed" validated nothing. These tests lock the fix at the level
that matters: a planted defect in an ASV run is CAUGHT, not skipped.

The banned-rank check (d) is the probe: inject an invented Linnean rank (Craniata as a Class) into
the ASV consensus and assert the ASV-run path FAILs on it — the exact defect the otu-blind verifier
could not see.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.verify import Report, _present_sets, check_banned_ranks, check_linnean_backbone

_ASV_TABLE = "ASV_ID\taptl64b_02\taptl64f_01\nASV1\t10\t20\nASV2\t30\t40\n"
_CONSENSUS_HEADER = "ASV_ID,Consensus_Phylum,Consensus_Class\n"


def _asv_run(tmp: Path, class_value: str) -> Path:
    (tmp / "analysis" / "dada2").mkdir(parents=True, exist_ok=True)
    (tmp / "analysis" / "dada2" / "asv_table.tsv").write_text(_ASV_TABLE, encoding="utf-8")
    db = tmp / "databases" / "asv"
    db.mkdir(parents=True, exist_ok=True)
    (db / "taxonomy_consensus.csv").write_text(
        _CONSENSUS_HEADER + f"ASV1,Chordata,{class_value}\nASV2,Chordata,Actinopteri\n", encoding="utf-8"
    )
    return tmp


class VerifyAsvTests(unittest.TestCase):
    def test_present_sets_detects_asv_from_its_abundance_table(self):
        with tempfile.TemporaryDirectory() as d:
            run = _asv_run(Path(d), "Actinopteri")
            self.assertEqual(_present_sets(run), ["asv"])

    def test_banned_rank_in_the_asv_consensus_is_caught(self):
        # Craniata is a banned invented Class — the defect the otu-blind verifier skipped.
        with tempfile.TemporaryDirectory() as d:
            run = _asv_run(Path(d), "Craniata")
            report = Report(quiet=True)
            check_banned_ranks(run, report, "asv")
            self.assertEqual(report.failed, 1, "asv-aware verify must FAIL on a banned rank in databases/asv/")
            status, name, detail = report.rows[0]
            self.assertEqual(status, "FAIL")
            self.assertIn("(d/asv)", name)
            self.assertIn("craniata", detail.lower())

    def test_a_clean_asv_consensus_passes(self):
        with tempfile.TemporaryDirectory() as d:
            run = _asv_run(Path(d), "Actinopteri")
            report = Report(quiet=True)
            check_banned_ranks(run, report, "asv")
            self.assertEqual(report.failed, 0)
            self.assertEqual(report.rows[0][0], "PASS")


def _linnean_run(tmp: Path, db_sources, sources_cell: str, phylum: str = "") -> Path:
    """A run with just what check (f) reads: progress.json's requested sources and the consensus."""
    run_config = {} if db_sources is None else {"db_sources": db_sources}
    (tmp / "progress.json").write_text(json.dumps({"run_config": run_config}), encoding="utf-8")
    db = tmp / "databases" / "asv"
    db.mkdir(parents=True, exist_ok=True)
    (db / "taxonomy_consensus.csv").write_text(
        "ASV_ID,Consensus_Sources,Consensus_Phylum,Consensus_Class\n"
        f"ASV1,{sources_cell},{phylum},\n",
        encoding="utf-8",
    )
    return tmp


class LinneanBackboneProvenanceTests(unittest.TestCase):
    """Check (f) must tell a config choice from a death — read from what was REQUESTED."""

    def _run_f(self, **kwargs):
        with tempfile.TemporaryDirectory() as d:
            run = _linnean_run(Path(d), **kwargs)
            report = Report(quiet=True)
            check_linnean_backbone(run, report, "asv")
            return report.rows[0]

    def test_no_linnean_source_requested_is_not_applicable(self):
        # SILVA+PR2-only: empty Phylum..Family is a marker choice, not a failure — must not alarm.
        status, name, detail = self._run_f(db_sources=["silva", "pr2"], sources_cell="SILVA|PR2")
        self.assertEqual(status, "SKIP")
        self.assertIn("(f/asv)", name)
        self.assertIn("config choice", detail)

    def test_requested_linnean_source_that_died_still_fails(self):
        # NCBI was asked for and never reached the consensus — the real case the alarm guards.
        status, _, detail = self._run_f(db_sources=["silva", "pr2", "ncbi"], sources_cell="SILVA|PR2")
        self.assertEqual(status, "FAIL")
        self.assertIn("NCBI", detail)

    def test_a_linnean_source_that_reached_the_consensus_passes(self):
        status, _, _ = self._run_f(db_sources=["ncbi"], sources_cell="NCBI", phylum="Chordata")
        self.assertEqual(status, "PASS")

    def test_unrecorded_provenance_falls_back_to_fail(self):
        # No db_sources recorded: cannot prove a Linnean source was NOT asked for → stay strict.
        status, _, _ = self._run_f(db_sources=None, sources_cell="SILVA")
        self.assertEqual(status, "FAIL")


if __name__ == "__main__":
    unittest.main()
