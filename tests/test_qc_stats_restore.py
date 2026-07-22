"""The QC KPI cards must survive a restored session — reconstructed from the golden qc tables.

`update_qc_stats` was store-only: it read the four KPIs (Total samples, Pass rate, Mean reads,
Mean GC) straight from `store-qc-outputs`. On a reloaded session that store is empty even though
the run is on disk, so the cards blanked to "—" — the same session-restore gap the other pages
closed, left here "by shape, not by axis" (the DEVLOG follow-up).

QC is the fifth golden layer, so the check is STRONG FORM: not merely "does not blank", but
"reconstructs the RIGHT values". The reconstruction is semantic (a third restore form): the store
shape isn't on disk, but the same per-sample numbers are — `qc_overview.csv` carries total_reads
and pct_gc per sample, and every sample present is one FastQC succeeded on (build_qc_metrics only
tabulates parseable output), so the run pass rate the store would show is total/total.

Sourced from the COMMITTED golden qc (golden/fixture/qc), so it runs on a clean clone.
"""

from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GOLDEN_QC = ROOT / "golden" / "fixture" / "qc" / "qc"  # the frozen qc/ subtree (metrics/*.csv)


class QcStatsRestoreTests(unittest.TestCase):
    def setUp(self):
        if not (GOLDEN_QC / "metrics" / "qc_overview.csv").exists():
            self.skipTest("qc golden not frozen")
        self._tmp = Path(tempfile.mkdtemp())
        shutil.copytree(GOLDEN_QC, self._tmp / "qc")

    def tearDown(self):
        shutil.rmtree(self._tmp, ignore_errors=True)

    def test_reconstructs_the_right_values_from_golden_qc(self):
        # Empty store = a restored session; run_config points at the on-disk (golden) run.
        from ui.callbacks_qc import _qc_stats

        stats = _qc_stats({}, {"output_dir": str(self._tmp)}, None)
        self.assertIsNotNone(stats, "restored session blanked instead of reconstructing from disk")
        self.assertEqual(stats["total"], 2)
        self.assertEqual(stats["pass_pct"], "100%")
        self.assertEqual(stats["mean_reads"], "5,000")
        self.assertEqual(stats["mean_gc"], "54%")
        self.assertTrue(stats["all_pass"])

    def test_reconstruction_matches_the_store_path_for_the_same_run(self):
        # The restore path must reproduce what the populated store would have shown, not a
        # near-miss: same total, pass rate, mean reads and mean GC.
        from ui.callbacks_qc import _qc_stats

        store = {
            "status": "success",
            "outputs": {
                "aptl64b_02": {"success": True, "metrics": {"total_sequences_r1": 5000, "pct_gc_r1": 54.0}},
                "aptl64f_01": {"success": True, "metrics": {"total_sequences_r1": 5000, "pct_gc_r1": 54.0}},
            },
        }
        from_store = _qc_stats(store, {}, None)
        from_disk = _qc_stats({}, {"output_dir": str(self._tmp)}, None)
        for key in ("total", "pass_pct", "mean_reads", "mean_gc"):
            self.assertEqual(from_store[key], from_disk[key], f"{key} diverged between store and restore")

    def test_no_store_and_no_disk_blanks(self):
        # With neither a store nor a resolvable run, the cards honestly show em-dashes.
        from ui.callbacks_qc import _qc_stats

        self.assertIsNone(_qc_stats({}, {"output_dir": str(self._tmp / "nope")}, None))


if __name__ == "__main__":
    unittest.main()
