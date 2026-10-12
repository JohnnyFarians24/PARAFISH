"""The full pipeline must resolve the run's marker the way Setup does: primer preset first.

The Setup page, the standalone Databases callback and ``RunConfig.marker`` all go through
``resolve_run_marker`` (preset first, file names as fallback). The full-pipeline thread —
the path behind "Start Pipeline" and ``scripts/run_pipeline.py`` — did not: it called
``resolve_marker`` over the samples' file-name markers alone. A COI/Leray run on files
named ``aptl64b_01`` therefore reached the Databases stage as ``OTHER``, and SILVA/PR2 —
SSU references that cannot carry COI — were not marked ``not_applicable``.

Two layers, on purpose. The helper tests pin the three resolution cases. The wiring test
drives the real ``_run_pipeline_thread`` and captures the ``marker`` it hands to the
taxonomy stage, because a correct helper that the thread does not call fixes nothing —
that test was seen to fail with ``'OTHER' != 'COI'`` before the fix.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from core.models import AnalysisBranch, AnalysisMode, Marker, RunConfig, Sample, SamplePaths, StepOutput
from databases.taxonomy_runner import _marker_incompatible


def _config(tmp: Path, marker: Marker, preset) -> RunConfig:
    r1 = tmp / "aptl64b_01_R1.fastq"
    r1.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    sample = Sample(
        internal_id="aptl64b_01",
        display_label="APTL64B",
        paths=SamplePaths(r1=r1),
        mode=AnalysisMode.SE,
        marker=marker,
    )
    config = RunConfig(run_id="t", samples=[sample], branch=AnalysisBranch.VSEARCH, output_dir=tmp / "run")
    if preset is not None:
        config.vsearch_params["primer_preset"] = preset
    return config


class PipelineRunMarkerTests(unittest.TestCase):
    def _marker(self, marker: Marker, preset) -> str:
        from ui.callbacks_progress import _pipeline_run_marker

        with tempfile.TemporaryDirectory() as tmp:
            return _pipeline_run_marker(_config(Path(tmp), marker, preset))

    def test_coi_preset_wins_over_unlabelled_file_names(self):
        marker = self._marker(Marker.OTHER, "COI_Leray")
        self.assertEqual(marker, "COI")
        # The consequence that matters: SSU references are refused before any work.
        self.assertTrue(_marker_incompatible("SILVA", marker))
        self.assertTrue(_marker_incompatible("PR2", marker))
        self.assertFalse(_marker_incompatible("BOLD", marker))

    def test_18s_files_with_the_auto_preset_still_resolve_to_18s(self):
        self.assertEqual(self._marker(Marker.SSU_18S, "auto"), "18S")
        self.assertEqual(self._marker(Marker.SSU_18S, None), "18S")

    def test_no_preset_and_no_marker_is_unchanged(self):
        from databases.reference_filter import resolve_marker

        for preset in ("", "auto", None):
            with self.subTest(preset=preset):
                # Exactly what the thread computed before: the file-name marker alone.
                self.assertEqual(self._marker(Marker.OTHER, preset), resolve_marker([Marker.OTHER]))


class PipelineThreadUsesTheRunMarkerTests(unittest.TestCase):
    def test_the_taxonomy_stage_receives_the_preset_marker(self):
        import ui.callbacks_progress as progress

        captured = {}

        def fake_taxonomy(fasta_path, output_dir, **kwargs):
            captured.update(kwargs)
            return StepOutput(step_name="taxonomy", success=False, files={}, error_message="stub")

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            config = _config(tmp_path, Marker.OTHER, "COI_Leray")
            config.enabled_stages = {"qc": False, "preprocess": False, "analysis": True, "databases": True, "phylo": False}
            otus = tmp_path / "otus.fasta"
            otus.write_text(">OTU1;size=3\nACGTACGT\n", encoding="utf-8")
            feature_set = {"key": "otu", "label": "OTU (VSEARCH)", "fasta": otus, "matrix": None}
            ok = StepOutput(step_name="stub", success=True, files={"filtered_fasta": otus, "otus_fasta": otus})
            # No real reference: this run is about the marker, and a run with the Databases
            # stage on fingerprints the references it finds — which, pointed at the
            # developer's .env, would hash the real SILVA and leave a cache file beside it.
            no_references = SimpleNamespace(silva=None, pr2=None, pr2_tax=None, pr2_dada2=None, bold_dir=None, blast=None)

            with patch.object(progress, "get_db_paths", return_value=no_references), \
                 patch.object(progress, "run_vsearch_filter_sample", return_value=ok), \
                 patch.object(progress, "run_vsearch_pool_and_cluster", return_value=ok), \
                 patch.object(progress, "feature_sets_from_outputs", return_value=[feature_set]), \
                 patch.object(progress, "build_metrics_from_outputs", return_value={}), \
                 patch.object(progress, "build_preprocess_metrics", return_value={}), \
                 patch.object(progress, "build_database_metrics", return_value={}), \
                 patch.object(progress, "_apply_cleanup_post_run", return_value=None), \
                 patch.object(progress, "run_taxonomy_pipeline", side_effect=fake_taxonomy):
                progress._run_pipeline_thread(config, {"db_sources": ["silva", "pr2"]}, config.output_dir / "progress.json")

        self.assertIn("marker", captured, "the taxonomy stage was never reached")
        self.assertEqual(captured["marker"], "COI")


if __name__ == "__main__":
    unittest.main()
