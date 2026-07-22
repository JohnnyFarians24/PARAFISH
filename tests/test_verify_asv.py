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

import tempfile
import unittest
from pathlib import Path

from scripts.verify import Report, _present_sets, check_banned_ranks

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


if __name__ == "__main__":
    unittest.main()
