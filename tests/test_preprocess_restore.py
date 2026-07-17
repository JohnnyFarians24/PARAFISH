"""A fresh session must reconstruct the Preprocess page from disk, not blank it.

The stage stores are per-session. Close the app and reopen tomorrow and store-run-config
is empty, so the three Preprocess render callbacks — which read ``output_dir`` straight
off that store — returned ``{}`` and showed "Run Preprocess...", even though the run is
still on disk. This is the second half of the Results-page bug (the disk path vs the
store path), and the OTU golden already pins the numbers the page must show, so this is
the STRONG form of the check: it asserts the page reconstructs AND the values match the
frozen golden, not merely "not blank".

Two independent gaps, each seen to fail before the fix:
 * LOGIC — ``_preprocess_metric_paths`` had no ``project_data`` and no ``latest_run_dir``
   fallback, so with empty stores it returned ``{}`` (it did not even accept a project).
 * WIRING — the three render callbacks did not declare ``store-project``, so the callback
   could not reach the run on disk even once the helper could.

And a trap the run-dir listing exposed: ``latest_run_dir`` picked ``max(..., key=mtime)``.
A test/junk run left in the project dir with a newer mtime but an older run_id would win,
so the restore would reconstruct from garbage and go green over it. The adversarial case
below pins that: newest *run_id* wins, not newest *mtime*.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FROZEN = ROOT / "golden" / "real" / "preprocess" / "preprocess" / "metrics" / "read_tracking.csv"
SRC_PREPROCESS = ROOT / "runs" / "_real" / "current" / "preprocess"


def _make_run(runs_root: Path, project_slug: str, run_id: str) -> Path:
    run_dir = runs_root / project_slug / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SRC_PREPROCESS, run_dir / "preprocess", dirs_exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, "done": True}), encoding="utf-8")
    return run_dir


class PreprocessRestoreLogicTests(unittest.TestCase):
    def setUp(self):
        self._prev = os.environ.get("RUNS_DIR")
        self.tmp = Path(tempfile.mkdtemp())
        os.environ["RUNS_DIR"] = str(self.tmp)

    def tearDown(self):
        if self._prev is None:
            os.environ.pop("RUNS_DIR", None)
        else:
            os.environ["RUNS_DIR"] = self._prev
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_reconstructs_frozen_numbers_from_disk_with_empty_session_stores(self):
        from ui.callbacks_preprocess import _preprocess_metric_paths

        _make_run(self.tmp, "restore_proj", "20260715T151700Z")
        # Fresh session/tab: no preprocess store, no run-config store, only store-project.
        paths = _preprocess_metric_paths(None, None, {"project_name": "restore_proj"})
        self.assertIn("read_tracking_csv", paths)
        self.assertEqual(Path(paths["read_tracking_csv"]).read_bytes(), FROZEN.read_bytes())

    def test_picks_projects_newest_run_id_not_a_newer_mtime_or_an_orphan(self):
        from ui.callbacks_preprocess import _preprocess_metric_paths

        # The correct run: newest run_id in the right project, real data.
        _make_run(self.tmp, "restore_proj", "20260715T151700Z")
        # A cross-project orphan with an even newer run_id — must be ignored (scoping).
        _make_run(self.tmp, "qc_validation_run", "20990101T000000Z")
        # An older-run_id decoy in the right project, GARBAGE data, forced to the newest
        # MTIME. mtime-based selection would pick this and reconstruct from garbage.
        decoy = _make_run(self.tmp, "restore_proj", "20200101T000000Z")
        (decoy / "preprocess" / "metrics" / "read_tracking.csv").write_text("sample\nGARBAGE\n", encoding="utf-8")
        future = os.stat(self.tmp).st_mtime + 100000
        os.utime(decoy, (future, future))

        paths = _preprocess_metric_paths(None, None, {"project_name": "restore_proj"})
        self.assertIn("read_tracking_csv", paths)
        self.assertEqual(Path(paths["read_tracking_csv"]).read_bytes(), FROZEN.read_bytes())


class PreprocessRestoreWiringTests(unittest.TestCase):
    def test_render_callbacks_depend_on_store_project(self):
        from app import create_app

        app = create_app()
        client = app.server.test_client()
        deps = client.get("/_dash-dependencies").get_json()
        targets = {
            "preprocess-tracking-table.children",
            "preprocess-stat-input.children",
            "preprocess-plots.children",
        }
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
