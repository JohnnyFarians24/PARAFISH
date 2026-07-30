"""The Analysis OTU stat card must survive a restored session — locked against regression.

Measured (2026-07-30): `_update_analysis_stats` is store-only — it counts OTUs from the FASTA
path carried in `store-analysis-outputs`. On a reloaded tab that store is empty, but
`restore_run_state` repopulates it from the run's progress.json, whose analysis record carries
the absolute FASTA path. So the card does NOT blank on restore — the same shape as
`update_qc_stats`, which "looked open and was closed". This is not a fix; it is the lock the
measurement earned: without it, a future change that stops `restore_run_state` handing back the
analysis record (or `_analysis_stat_values` reading it) blanks the card silently, and the suite
would not notice — exactly the regression the qc-restore test exists to stop, one page over.

The test drives the REAL restore chain (restore_run_state → `_analysis_stat_values`) against a
self-contained tmp run, so it runs on a clean clone (the golden commits otu_matrix.csv but not
the OTU FASTA, and runs/_real is gitignored, so the count is pinned by a fixture FASTA here).
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

_FASTA = ">OTU1;size=10\nACGT\n>OTU2;size=5\nTTAA\n>OTU3;size=2\nGGCC\n"  # 3 OTUs


class AnalysisStatsRestoreTests(unittest.TestCase):
    def setUp(self):
        self._prev = os.environ.get("RUNS_DIR")
        self.tmp = Path(tempfile.mkdtemp())
        os.environ["RUNS_DIR"] = str(self.tmp)
        self.project_data = {"project_name": "stats_restore_proj"}
        run_id = "20260730T101500Z"
        run_dir = self.tmp / "stats_restore_proj" / run_id
        (run_dir / "analysis" / "pooled").mkdir(parents=True, exist_ok=True)
        fasta = run_dir / "analysis" / "pooled" / "otus.fasta"
        fasta.write_text(_FASTA, encoding="utf-8")
        # A FINISHED run whose analysis record carries the absolute FASTA path — the shape
        # restore_run_state hands back to store-analysis-outputs on a fresh tab.
        (run_dir / "progress.json").write_text(
            json.dumps({
                "run_id": run_id,
                "done": True,
                "analysis_record": {
                    "status": "success",
                    "outputs": {"combined_fasta_vsearch": str(fasta)},
                },
            }),
            encoding="utf-8",
        )

    def tearDown(self):
        if self._prev is None:
            os.environ.pop("RUNS_DIR", None)
        else:
            os.environ["RUNS_DIR"] = self._prev
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_restored_session_reports_the_otu_count_not_a_dash(self):
        from ui.callbacks_analysis import _analysis_stat_values
        from ui.runtime_state import restore_run_state

        restored = restore_run_state(self.project_data)
        self.assertIsNotNone(restored, "restore_run_state found nothing to restore")
        values = _analysis_stat_values(restored["analysis"], restored["run_config"])
        self.assertEqual(values["otu"], 3, "restored session blanked the OTU count instead of reading the run's FASTA")

    def test_empty_store_and_no_run_blanks_honestly(self):
        # With neither a store nor a resolvable FASTA path, the card shows an em-dash, not a zero.
        from ui.callbacks_analysis import _analysis_stat_values

        self.assertEqual(_analysis_stat_values({}, {})["otu"], "—")


if __name__ == "__main__":
    unittest.main()
