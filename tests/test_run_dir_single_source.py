"""Task 7 — one selector for "the newest run to show", shared by Results and Eixo B.

Two functions named ``latest_run_dir`` answered the same question independently:
``ui.runtime_state.latest_run_dir`` (Eixo B page fallback) was hardened to order by
run_id and skip incomplete runs; ``core.run_registry.latest_run_dir`` (Results, via
``run_dir_for``) still ordered by mtime with no completeness gate. Same responsibility,
one hardened and one not — the phylo-guard class.

It has a concrete consequence, not just debt: with an OOM orphan that is the newest run,
Results (mtime, no gate) picks the orphan while the Eixo B pages (run_id, complete-only)
pick the last finished run — two different runs in ONE session, worse than either being
wrong alone. The fix collapses the SELECTION onto one primitive
(``core.run_registry.newest_complete_run``); the legitimately different project scoping
(Results falls back across projects, Eixo B is project-scoped) is kept.

Seen to fail before the fix: with run_registry still on mtime + no gate, Results and Eixo
B choose different runs and these assertions go red.
"""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


def _write_run(root: Path, project: str, run_id: str, state: dict, mtime: float | None = None) -> Path:
    run_dir = root / project / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    progress = run_dir / "progress.json"
    progress.write_text(json.dumps({"run_id": run_id, **state}), encoding="utf-8")
    if mtime is not None:
        os.utime(progress, (mtime, mtime))
        os.utime(run_dir, (mtime, mtime))
    return run_dir


_COMPLETE = {"done": True, "stages": {"QC": "success"}}
_ORPHAN = {"stages": {"QC": "success", "Analysis": "running"}}

# active_id/projects shape (what run_registry reads) whose name is already slug-like so
# both selectors resolve to the same on-disk project directory.
_PROJECT = {"active_id": "p", "projects": {"p": {"project_name": "proj"}}}


class ResultsAndEixoBAgreeTests(unittest.TestCase):
    def setUp(self):
        self._prev = os.environ.get("RUNS_DIR")
        self.tmp = Path(tempfile.mkdtemp())
        os.environ["RUNS_DIR"] = str(self.tmp)

    def tearDown(self):
        if self._prev is None:
            os.environ.pop("RUNS_DIR", None)
        else:
            os.environ["RUNS_DIR"] = self._prev
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def _results_pick(self):
        from core.run_registry import run_dir_for

        # No session output_dir -> exercises the disk fallback (the shared question).
        with patch("core.run_registry.get_runs_dir", lambda: self.tmp):
            return run_dir_for(None, _PROJECT)

    def _eixo_b_pick(self):
        from ui.runtime_state import latest_run_dir

        return latest_run_dir(_PROJECT)

    def test_they_agree_when_the_newest_run_is_an_orphan(self):
        """The false-divergence guard. The newest run (by run_id AND mtime) is an OOM
        orphan; both surfaces must fall to the older run that finished."""
        complete = _write_run(self.tmp, "proj", "20260101T000000Z", _COMPLETE, mtime=1_000)
        _write_run(self.tmp, "proj", "20260102T000000Z", _ORPHAN, mtime=2_000)  # newest, orphan

        results = self._results_pick()
        eixo_b = self._eixo_b_pick()
        self.assertEqual(results, eixo_b, "Results and Eixo B must pick the same run")
        self.assertEqual(eixo_b, complete, "both must pick the finished run, not the orphan")

    def test_ordering_is_by_run_id_not_mtime(self):
        """An older run_id with a freshly bumped mtime (a re-freeze, a copy) must not win
        over the genuinely newest run_id."""
        _write_run(self.tmp, "proj", "20260101T000000Z", _COMPLETE, mtime=9_000)  # bumped mtime
        newest = _write_run(self.tmp, "proj", "20260103T000000Z", _COMPLETE, mtime=1_000)

        self.assertEqual(self._results_pick(), newest, "Results must order by run_id, not mtime")
        self.assertEqual(self._eixo_b_pick(), newest)


class SharedPrimitiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_newest_complete_run_skips_orphans_and_orders_by_run_id(self):
        from core.run_registry import newest_complete_run

        _write_run(self.tmp, "proj", "20260102T000000Z", _ORPHAN, mtime=5_000)
        complete = _write_run(self.tmp, "proj", "20260101T000000Z", _COMPLETE, mtime=1_000)
        run_dirs = [p.parent for p in (self.tmp / "proj").glob("*/progress.json")]

        self.assertEqual(newest_complete_run(run_dirs), complete)

    def test_newest_complete_run_returns_none_when_all_orphans(self):
        from core.run_registry import newest_complete_run

        _write_run(self.tmp, "proj", "20260101T000000Z", _ORPHAN)
        run_dirs = [p.parent for p in (self.tmp / "proj").glob("*/progress.json")]
        self.assertIsNone(newest_complete_run(run_dirs))

    def test_both_latest_run_dir_use_the_same_primitive(self):
        """Structural: neither selector may re-derive completeness/ordering on its own."""
        registry_src = (Path(__file__).resolve().parent.parent / "core" / "run_registry.py").read_text()
        runtime_src = (Path(__file__).resolve().parent.parent / "ui" / "runtime_state.py").read_text()
        self.assertIn("def newest_complete_run", registry_src)
        self.assertIn("newest_complete_run", runtime_src)
        # The old mtime ordering must be gone from the run selector.
        self.assertNotIn("st_mtime", registry_src, "run selection must not order by mtime")


if __name__ == "__main__":
    unittest.main()
