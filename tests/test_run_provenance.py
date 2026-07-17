"""The config a run was made under must be recorded in its own progress.json.

progress.json used to carry only per-stage status, so two runs that differ ONLY in
configuration — which reference DBs were queried, ``min_unique_size``, whether primers
were stripped — produced different outputs that read as a regression against golden
when they were nothing of the sort. This check pins the fix: two fixture runs that
differ only in ``--db-sources`` must be distinguishable by the provenance recorded in
their *real* progress.json — the artifact the pipeline writes, not a hand-built one.

Anchored on the real pipeline entrypoint: it drives ``build_stores`` (the same helper
``scripts/run_pipeline.py`` uses) and the real ``_run_pipeline_thread``. Every stage is
disabled so no external tool is required — the run_config block is written at init,
before any stage, which is exactly why it survives a run whose tools never ran.

Seen to fail before the fix: progress.json had no ``run_config`` key, so both runs were
indistinguishable and the ``progress["run_config"]`` lookup raised ``KeyError``.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _load_run_pipeline():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    spec = importlib.util.spec_from_file_location("run_pipeline", ROOT / "scripts" / "run_pipeline.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_and_read_provenance(tmp: Path, db_sources: list) -> dict:
    run_pipeline = _load_run_pipeline()
    from ui.callbacks_progress import _run_pipeline_thread
    from ui.runtime_state import run_config_from_store

    out_dir = tmp / ("_".join(db_sources) or "none")
    out_dir.mkdir(parents=True, exist_ok=True)
    # stages=[] -> every stage disabled: no external tool runs, but the run_config
    # block is written at init, so this exercises the real pipeline offline.
    project_data, samples_data, run_config_data = run_pipeline.build_stores(
        run_pipeline.FIXTURE_INPUTS, out_dir, [], db_sources, "testrun", {"threads": 1}
    )
    run_config = run_config_from_store(project_data, samples_data, run_config_data)
    run_config.run_id = "testrun"
    run_config.output_dir = out_dir
    _run_pipeline_thread(run_config, run_config_data, out_dir / "progress.json")
    return json.loads((out_dir / "progress.json").read_text())["run_config"]


class RunProvenanceTests(unittest.TestCase):
    def test_db_sources_are_recorded_and_distinguish_runs(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            prov_two = _run_and_read_provenance(tmp, ["silva", "pr2"])
            prov_one = _run_and_read_provenance(tmp, ["silva"])
        self.assertEqual(prov_two["db_sources"], ["silva", "pr2"])
        self.assertEqual(prov_one["db_sources"], ["silva"])
        self.assertNotEqual(prov_two["db_sources"], prov_one["db_sources"])

    def test_provenance_captures_the_other_config_axes(self):
        with tempfile.TemporaryDirectory() as d:
            prov = _run_and_read_provenance(Path(d), ["silva"])
        # The axes a golden diff must be read under, all resolvable before any stage.
        self.assertEqual(prov["markers"], ["18S"])
        self.assertIn("min_unique_size", prov)
        self.assertIn("primers", prov)
        self.assertIn(prov["ncbi_mode"], ("local", "remote"))


if __name__ == "__main__":
    unittest.main()
