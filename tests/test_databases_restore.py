"""A fresh session must reconstruct the Databases page from disk, not blank it.

Same class as the Preprocess page (Eixo B): the stage stores are per-session, so on a
restart or a fresh tab ``store-run-config`` is empty and the render callbacks, which read
``output_dir`` straight off it, returned ``{}`` / "Run the Databases stage...". The
taxonomy golden pins the numbers, so this is the STRONG form — reconstructs AND the bytes
match golden — on the SELECTOR-FREE callbacks only.

Deliberately clear of the feature-set selector (Eixo A). The strong-form path resolves
the ``otu`` tables through ``_database_metric_paths`` / ``_taxonomy_file`` with the plain
``otu`` default — the ``databases-feature-set`` widget is nowhere in it. Proving those two
selector callbacks (``_render_hits``, ``_render_consensus``) reconstruct the RIGHT set is
the "otu-hardcoded" Eixo A problem, unprovable without an ASV golden, and is left out.

Seen to fail before the fix:
 * ``_database_metric_paths`` / ``_resolve_output_dir`` had no ``project_data`` and no
   ``latest_run_dir`` fallback → ``{}``.
 * the selector-free render callbacks did not declare ``store-project``.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FROZEN_RESOLUTION = ROOT / "golden" / "real" / "taxonomy" / "databases" / "otu" / "metrics" / "taxonomic_resolution.csv"
FROZEN_CONFLICT = ROOT / "golden" / "real" / "taxonomy" / "databases" / "otu" / "conflict_report.csv"
SRC_DATABASES = ROOT / "runs" / "_real" / "current" / "databases"


def _make_run(runs_root: Path, project_slug: str, run_id: str) -> Path:
    run_dir = runs_root / project_slug / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SRC_DATABASES, run_dir / "databases", dirs_exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, "done": True}), encoding="utf-8")
    return run_dir


class DatabasesRestoreLogicTests(unittest.TestCase):
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

    def test_metric_paths_reconstructs_frozen_resolution_with_empty_stores(self):
        from ui.callbacks_databases import _database_metric_paths

        paths = _database_metric_paths(None, None, self.project_data)
        self.assertIn("resolution_csv", paths)
        self.assertEqual(Path(paths["resolution_csv"]).read_bytes(), FROZEN_RESOLUTION.read_bytes())

    def test_resolve_output_dir_feeds_taxonomy_file_to_frozen_conflict_report(self):
        from ui.callbacks_databases import _resolve_output_dir, _taxonomy_file

        run_output_dir = _resolve_output_dir(None, self.project_data)
        self.assertIsNotNone(run_output_dir)
        path = _taxonomy_file({}, run_output_dir, "conflict_report_csv", "conflict_report.csv")
        self.assertIsNotNone(path)
        self.assertEqual(Path(path).read_bytes(), FROZEN_CONFLICT.read_bytes())


class DatabasesRestoreWiringTests(unittest.TestCase):
    def test_selector_free_render_callbacks_depend_on_store_project(self):
        from app import create_app

        app = create_app()
        client = app.server.test_client()
        deps = client.get("/_dash-dependencies").get_json()
        targets = {
            "databases-resolution.children",
            "databases-confidence.children",
            "databases-conflicts.children",
            "databases-contaminants.children",
            "databases-prevalence.children",
            "databases-feature-comparison.children",
            "databases-visuals.children",
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
