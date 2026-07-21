"""_render_hits / _render_consensus: the two SELECTOR callbacks deferred until an ASV golden.

``test_databases_restore.py`` restored the Databases page's selector-FREE callbacks in the
strong form (reconstruct + bytes match golden), but left ``_render_hits`` and
``_render_consensus`` out on purpose: proving a selector callback reconstructs the RIGHT set is
the Eixo A "otu-hardcoded" problem, and it was unprovable without an ASV taxonomy golden. That
golden now exists (``golden/fixture_asv/taxonomy``, SINTAX-seeded, determinism-gated), so the
proof is finally possible and this is it.

Two halves, same as occupancy:
 * Eixo B restore — both callbacks read ``output_dir`` straight off the per-session
   ``store-run-config``, so a fresh tab blanked them; they now resolve through
   ``_resolve_output_dir`` with a ``store-project`` fallback.
 * Eixo A right-set — with the run reconstructed from disk, the ``asv`` selector yields the asv
   tables byte-identical to golden, and the ``otu`` selector on an asv-only run yields NOTHING
   (the set is honored, not defaulted to otu).

Seen to fail before the fix: the two render callbacks did not declare ``store-project``.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FROZEN_ASV_HITS = ROOT / "golden" / "fixture_asv" / "taxonomy" / "databases" / "asv" / "taxonomy_hits.csv"
FROZEN_ASV_CONSENSUS = ROOT / "golden" / "fixture_asv" / "taxonomy" / "databases" / "asv" / "taxonomy_consensus.csv"
SRC_DATABASES = ROOT / "runs" / "_fixture_asv" / "current" / "databases"


def _make_run(runs_root: Path, project_slug: str, run_id: str) -> Path:
    run_dir = runs_root / project_slug / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SRC_DATABASES, run_dir / "databases", dirs_exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, "done": True}), encoding="utf-8")
    return run_dir


class HitsConsensusRestoreLogicTests(unittest.TestCase):
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

    def test_reconstructed_dir_feeds_selected_sets_hits(self):
        from ui.callbacks_databases import _load_hits_frame, _resolve_output_dir
        from ui.feature_sets import taxonomy_file

        run_output_dir = _resolve_output_dir(None, self.project_data)
        self.assertIsNotNone(run_output_dir)

        # selector = asv → the asv hits table, byte-identical to golden, and the loader reads it.
        asv_path = taxonomy_file({}, run_output_dir, "asv", "taxonomy_hits_csv")
        self.assertIsNotNone(asv_path)
        self.assertEqual(Path(asv_path).read_bytes(), FROZEN_ASV_HITS.read_bytes())
        frame, error = _load_hits_frame({}, run_output_dir, "asv", 0, 0, 0)
        self.assertIsNone(error)
        self.assertIsNotNone(frame)
        self.assertFalse(frame.empty)

        # selector = otu on an asv-only run → nothing (the set is honored, not defaulted).
        self.assertIsNone(taxonomy_file({}, run_output_dir, "otu", "taxonomy_hits_csv"))
        otu_frame, otu_error = _load_hits_frame({}, run_output_dir, "otu", 0, 0, 0)
        self.assertIsNone(otu_frame)
        self.assertIsNone(otu_error)

    def test_reconstructed_dir_feeds_selected_sets_consensus(self):
        from ui.callbacks_databases import _find_consensus_csv, _resolve_output_dir

        run_output_dir = _resolve_output_dir(None, self.project_data)
        self.assertIsNotNone(run_output_dir)

        # selector = asv → the asv consensus table, byte-identical to golden.
        asv_path = _find_consensus_csv({}, run_output_dir, "asv")
        self.assertIsNotNone(asv_path)
        self.assertEqual(Path(asv_path).read_bytes(), FROZEN_ASV_CONSENSUS.read_bytes())

        # selector = otu on an asv-only run → nothing (the set is honored, not defaulted).
        self.assertIsNone(_find_consensus_csv({}, run_output_dir, "otu"))


class HitsConsensusWiringTests(unittest.TestCase):
    def test_selector_render_callbacks_depend_on_store_project(self):
        from app import create_app

        app = create_app()
        client = app.server.test_client()
        deps = client.get("/_dash-dependencies").get_json()
        targets = {"databases-hits.children", "databases-consensus.children"}
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
