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
        # The consensus is a DATABASES-stage output, in its per-feature-set (taxonomy_by_set)
        # location on disk, carried by db_record — not by the analysis record.
        (run_dir / "databases" / "otu").mkdir(parents=True, exist_ok=True)
        consensus = run_dir / "databases" / "otu" / "taxonomy_consensus.csv"
        consensus.write_text("feature_id,Consensus_Genus\nOTU1,Diplostomum\n", encoding="utf-8")
        # A FINISHED run whose stage records carry the absolute paths — the shape
        # restore_run_state hands back to the session stores on a fresh tab.
        (run_dir / "progress.json").write_text(
            json.dumps({
                "run_id": run_id,
                "done": True,
                "analysis_record": {
                    "status": "success",
                    "outputs": {"combined_fasta_vsearch": str(fasta)},
                },
                "db_record": {
                    "status": "success",
                    "outputs": {"taxonomy_by_set": {"otu": {"files": {"taxonomy_consensus_csv": str(consensus)}}}},
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

    def test_restored_session_lights_the_consensus_card_from_the_databases_store(self):
        # The crossed wire: before the fix the card read store-ANALYSIS-outputs, which never
        # carries the consensus, so this branch could never fire even with the file on disk
        # and db_record repopulated. Seen to fail against that code, "—" != "ready".
        from ui.callbacks_analysis import _analysis_stat_values
        from ui.runtime_state import restore_run_state

        restored = restore_run_state(self.project_data)
        values = _analysis_stat_values(
            restored["analysis"], restored["run_config"], restored["databases"], self.project_data
        )
        self.assertEqual(values["consensus"], "ready", "consensus exists + db_record repopulated, yet the card stayed pending")

    def test_flat_pre_feature_set_consensus_shape_also_lights(self):
        # A run from before feature sets kept a single FLAT taxonomy.files entry (no
        # taxonomy_by_set). The card must light for that shape too, via the same canonical
        # resolver — reading only one shape would light for OTU and stay dead elsewhere.
        from ui.callbacks_analysis import _analysis_stat_values

        flat = self.tmp / "flat_consensus.csv"
        flat.write_text("feature_id\nOTU1\n", encoding="utf-8")
        db_store = {"outputs": {"taxonomy": {"files": {"taxonomy_consensus_csv": str(flat)}}}}
        self.assertEqual(_analysis_stat_values({}, {}, db_store, None)["consensus"], "ready")

    def test_consensus_stays_pending_when_none_exists(self):
        from ui.callbacks_analysis import _analysis_stat_values

        values = _analysis_stat_values({}, {"output_dir": str(self.tmp / "nope")}, {}, None)
        self.assertEqual(values["consensus"], "—")


if __name__ == "__main__":
    unittest.main()
