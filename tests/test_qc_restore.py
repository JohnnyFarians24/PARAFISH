"""A fresh session must reconstruct the QC page from disk, not blank it.

Eixo B, page 4/4 — the last, and the cleanest: QC is pre-clustering and per-sample, so it
has NO feature-set selector and no Eixo A entanglement at all. The whole disk-reconstructing
surface is provable now that the QC golden exists (frozen this session).

Two render callbacks reconstruct from disk via ``_qc_metric_paths`` and blank on a fresh
tab: ``_render_qc_plots`` (reads ``qc_length_distribution.csv`` / ``qc_per_base_quality.csv``)
and ``_render_qc_sample_table`` (reads ``qc_overview.csv``). ``update_qc_stats`` is store-only
by design (it reads the per-sample outputs dict, not a metrics CSV) — a distinct restore
shape, left out and noted, not in scope here.

Seen to fail before the fix: ``_qc_metric_paths`` had no ``project_data`` / ``latest_run_dir``
fallback (→ ``{}``), and the two render callbacks did not declare ``store-project``.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FROZEN_LENGTH = ROOT / "golden" / "real" / "qc" / "qc" / "metrics" / "qc_length_distribution.csv"
FROZEN_OVERVIEW = ROOT / "golden" / "real" / "qc" / "qc" / "metrics" / "qc_overview.csv"
SRC_QC = ROOT / "runs" / "_real" / "current" / "qc"


def _make_run(runs_root: Path, project_slug: str, run_id: str) -> Path:
    run_dir = runs_root / project_slug / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SRC_QC, run_dir / "qc", dirs_exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, "done": True}), encoding="utf-8")
    return run_dir


class QCRestoreLogicTests(unittest.TestCase):
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

    def test_reconstructs_frozen_qc_tables_with_empty_session_stores(self):
        from ui.callbacks_qc import _qc_metric_paths

        paths = _qc_metric_paths(None, None, self.project_data)
        self.assertIn("qc_length_distribution_csv", paths)
        self.assertIn("qc_overview_csv", paths)
        self.assertEqual(Path(paths["qc_length_distribution_csv"]).read_bytes(), FROZEN_LENGTH.read_bytes())
        self.assertEqual(Path(paths["qc_overview_csv"]).read_bytes(), FROZEN_OVERVIEW.read_bytes())


class QCRestoreWiringTests(unittest.TestCase):
    def test_render_callbacks_depend_on_store_project(self):
        from app import create_app

        app = create_app()
        client = app.server.test_client()
        deps = client.get("/_dash-dependencies").get_json()
        targets = {"qc-plots.children", "qc-sample-table.children"}
        checked = set()
        for cb in deps:
            out = cb["output"]
            outs = [o.strip(".") for o in out.strip(".").split("...")] if out.startswith("..") else [out]
            for target in targets:
                if target in outs:
                    ids = {d["id"] for d in cb["inputs"]} | {d["id"] for d in cb["state"]}
                    self.assertIn("store-project", ids, f"{target}: render callback cannot see store-project")
                    checked.add(target)
        self.assertEqual(checked, targets, f"target callbacks not found: {targets - checked}")


if __name__ == "__main__":
    unittest.main()
