"""_render_occupancy: the deferred-whole callback whose scope, measured, was only restore.

The registry filed occupancy as deferred-whole ("Eixo B restore + right-set proof"). Measured
against the fixture_asv anchor, the right-set half was never broken: occupancy already passes
its ``analysis-occupancy-set`` selector value straight to ``analysis_metrics_paths(dir, set)``.
The only gap is the Eixo B restore — it read ``output_dir`` off the session store with no
``store-project`` fallback, so a fresh tab blanked it. So this is a plain restore wiring, and
the proof is: the reconstructed run dir feeds the selected set's occupancy table, byte-for-byte
against ``golden/fixture_asv`` for the asv path.

Wiring seen to fail before the fix: the callback did not declare ``store-project``.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FROZEN_ASV_OCCUPANCY = (
    ROOT / "golden" / "fixture_asv" / "analysis" / "analysis" / "metrics" / "asv" / "feature_occupancy.csv"
)
SRC_ANALYSIS = ROOT / "runs" / "_fixture_asv" / "current" / "analysis"


def _make_run(runs_root: Path, project_slug: str, run_id: str) -> Path:
    run_dir = runs_root / project_slug / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SRC_ANALYSIS, run_dir / "analysis", dirs_exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, "done": True}), encoding="utf-8")
    return run_dir


class OccupancyRestoreLogicTests(unittest.TestCase):
    def setUp(self):
        self._prev = os.environ.get("RUNS_DIR")
        self.tmp = Path(tempfile.mkdtemp())
        os.environ["RUNS_DIR"] = str(self.tmp)
        _make_run(self.tmp, "restore_proj", "20260715T151700Z")
        self.project_data = {"project_name": "restore_proj"}

    def tearDown(self):
        if self._prev is None:
            os.environ.pop("RUNS_DIR", None)
        else:
            os.environ["RUNS_DIR"] = self._prev
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_reconstructed_dir_feeds_selected_sets_occupancy(self):
        from analysis.analysis_metrics import analysis_metrics_paths
        from ui.callbacks_analysis import _resolve_output_dir

        run_output_dir = _resolve_output_dir(None, self.project_data)
        self.assertIsNotNone(run_output_dir)
        # selector = asv → the asv occupancy table, byte-identical to golden.
        asv = analysis_metrics_paths(Path(run_output_dir) / "analysis", "asv").get("occupancy_csv")
        self.assertIsNotNone(asv)
        self.assertEqual(Path(asv).read_bytes(), FROZEN_ASV_OCCUPANCY.read_bytes())
        # selector = otu on an asv-only run → nothing (the set is honored, not defaulted).
        self.assertIsNone(analysis_metrics_paths(Path(run_output_dir) / "analysis", "otu").get("occupancy_csv"))


class OccupancyWiringTests(unittest.TestCase):
    def test_occupancy_callback_depends_on_store_project(self):
        from app import create_app

        app = create_app()
        client = app.server.test_client()
        deps = client.get("/_dash-dependencies").get_json()
        for cb in deps:
            out = cb["output"]
            outs = [o.strip(".") for o in out.strip(".").split("...")] if out.startswith("..") else [out]
            if "analysis-occupancy.children" in outs:
                ids = {d["id"] for d in cb["inputs"]} | {d["id"] for d in cb["state"]}
                self.assertIn("store-project", ids, "occupancy cannot see store-project")
                return
        self.fail("analysis-occupancy callback not found")


if __name__ == "__main__":
    unittest.main()
