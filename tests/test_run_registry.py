"""The run registry: the UI's facts come from disk, not from session state."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from core.run_registry import (
    PER_SAMPLE_STAGES,
    list_run_samples,
    load_run_state,
    overall_status,
    run_dir_for,
    run_stages,
)


def _write_run(root: Path, *, analysis_status: str = "success") -> Path:
    run_dir = root / "proj" / "20260714T120000Z"
    run_dir.mkdir(parents=True)
    (run_dir / "progress.json").write_text(
        json.dumps(
            {
                "run_id": "20260714T120000Z",
                "done": True,
                "qc_record": {"status": "success", "outputs": {"s1": {"files": {}}, "s2": {"files": {}}}},
                "preprocess_record": {
                    "status": "success",
                    "outputs": {"s1": {"trimmomatic": {"files": {}}}, "s2": {"trimmomatic": {"files": {}}}},
                    "next_samples": [
                        {"internal_id": "s1", "display_label": "APTL64B"},
                        {"internal_id": "s2", "display_label": "APTL64F"},
                    ],
                },
                "analysis_record": {
                    "status": analysis_status,
                    "outputs": {
                        "s1": {"vsearch_filter": {"files": {}}},
                        # Stage-level artifacts sit in the same dict; they are not samples.
                        "combined_fasta": "/tmp/otus.fasta",
                        "feature_sets": [{"key": "otu"}],
                        "vsearch_pool": {"files": {}},
                    },
                },
                "db_record": {"status": "failure", "outputs": {}},
                "phylo_record": {"status": "skipped", "outputs": {}},
            }
        ),
        encoding="utf-8",
    )
    return run_dir


class DiskIsTheSourceOfTruthTests(unittest.TestCase):
    def test_the_run_is_found_with_no_ui_state_at_all(self):
        """The exact situation that killed the selectors: empty session stores.

        The dropdowns were driven by ``storage_type="session"`` stores that only the
        progress poller ever wrote, and the poller stops when the run ends — so a
        restart, a new tab, or simply arriving at Results later left them empty.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            _write_run(root)
            with patch("core.run_registry.get_runs_dir", lambda: root):
                run_dir = run_dir_for(None, None)  # no run-config, no project: nothing
                self.assertIsNotNone(run_dir)
                self.assertTrue(load_run_state(run_dir))

    def test_no_runs_on_disk_returns_none_rather_than_guessing(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch("core.run_registry.get_runs_dir", lambda: Path(temp_dir)):
                self.assertIsNone(run_dir_for(None, None))


class StageListingTests(unittest.TestCase):
    def test_failed_and_skipped_stages_are_still_listed(self):
        """Listing only successful stages is why the dropdown was empty on the run the
        user actually needed to debug: analysis and databases had failed."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            run_dir = _write_run(root, analysis_status="failure")
            stages = run_stages(load_run_state(run_dir))

        by_key = {stage["key"]: stage for stage in stages}
        self.assertEqual(len(stages), 5)  # every stage, always
        self.assertEqual(by_key["analysis"]["status"], "failure")
        self.assertEqual(by_key["databases"]["status"], "failure")
        self.assertEqual(by_key["phylo"]["status"], "skipped")

    def test_pooled_stages_are_marked_as_having_no_per_sample_output(self):
        self.assertIn("qc", PER_SAMPLE_STAGES)
        self.assertIn("preprocess", PER_SAMPLE_STAGES)
        # Databases and Phylo classify/build on the pooled feature set.
        self.assertNotIn("databases", PER_SAMPLE_STAGES)
        self.assertNotIn("phylo", PER_SAMPLE_STAGES)


class SampleListingTests(unittest.TestCase):
    def test_samples_come_with_their_display_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            run_dir = _write_run(Path(temp_dir))
            samples = list_run_samples(load_run_state(run_dir), "preprocess")

        self.assertEqual(samples, [{"label": "APTL64B", "value": "s1"}, {"label": "APTL64F", "value": "s2"}])

    def test_stage_level_artifacts_are_not_mistaken_for_samples(self):
        """``combined_fasta``, ``feature_sets`` and friends live in the same dict as the
        per-sample records, and naively iterating it offered them as selectable samples."""
        with tempfile.TemporaryDirectory() as temp_dir:
            run_dir = _write_run(Path(temp_dir))
            samples = list_run_samples(load_run_state(run_dir), "analysis")

        values = {entry["value"] for entry in samples}
        self.assertEqual(values, {"s1"})
        for artifact in ("combined_fasta", "feature_sets", "vsearch_pool"):
            self.assertNotIn(artifact, values)


class OverallStatusTests(unittest.TestCase):
    def test_a_half_good_run_is_partial_not_failed(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            run_dir = _write_run(Path(temp_dir))  # qc/pre/analysis ok, databases failed
            self.assertEqual(overall_status(load_run_state(run_dir)), "partial")

    def test_no_state_is_pending(self):
        self.assertEqual(overall_status({}), "pending")


if __name__ == "__main__":
    unittest.main()
