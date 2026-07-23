"""Task 6 — orphan runs must not masquerade as the active run and blank the restore.

A pipeline runs in a thread inside the server process, so a run left "running" on disk
is an OOM-killed / abandoned orphan (found: 7 such runs in runs/). Two protections, and
the second is the one that guards the restore:

  1. Startup reconciliation marks a "running but no live thread" run as interrupted, so
     it stops lying "running" forever.
  2. latest_run_dir — shared by Project and the four Eixo B pages — returns only a run
     that actually finished. If the NEWEST run is an orphan, restoring the Eixo B pages
     from it renders them BLANK, indistinguishable from the Eixo B restore bug having
     regressed. So latest_run_dir skips the orphan for the newest COMPLETE run, and when
     a project has no complete run it returns None (honest empty-state, not a blank
     orphan).

Decision for the all-orphans case: return None. A None makes the Eixo B callbacks show
their empty-state (they already handle it); returning the orphan would hand restore a run
whose tables never finished, i.e. the exact blank-that-looks-broken this task removes.

Seen to fail before the fix: with latest_run_dir picking the newest run_id regardless of
completion, test_prefers_the_newest_complete_run_over_a_newer_orphan gets the orphan and
goes red for the documented reason.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path


def _write_run(root: Path, project: str, run_id: str, state: dict) -> Path:
    run_dir = root / project / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, **state}), encoding="utf-8")
    return run_dir


_COMPLETE = {"done": True, "stages": {"QC": "success", "Analysis": "success"}}
_ORPHAN = {"stages": {"QC": "success", "Analysis": "running"}}


class ReconciliationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _state(self, run_dir: Path) -> dict:
        return json.loads((run_dir / "progress.json").read_text(encoding="utf-8"))

    def test_marks_a_running_orphan_as_interrupted(self):
        from core.run_registry import reconcile_interrupted_runs

        run = _write_run(self.tmp, "proj", "20260101T000000Z", _ORPHAN)
        touched = reconcile_interrupted_runs(runs_root=self.tmp)

        self.assertEqual(touched, [run])
        state = self._state(run)
        self.assertTrue(state.get("interrupted"))
        self.assertEqual(state["overall"], "interrupted")
        self.assertEqual(state["stages"]["Analysis"], "interrupted", "the running badge must stop lying")
        self.assertEqual(state["stages"]["QC"], "success", "terminal stages are left as they were")
        self.assertFalse(state.get("done"), "an interrupted run did NOT finish; done must stay unset")

    def test_leaves_a_finished_run_untouched(self):
        from core.run_registry import reconcile_interrupted_runs

        run = _write_run(self.tmp, "proj", "20260101T000000Z", _COMPLETE)
        before = self._state(run)
        touched = reconcile_interrupted_runs(runs_root=self.tmp)

        self.assertEqual(touched, [])
        self.assertEqual(self._state(run), before, "a completed run must not be rewritten")

    def test_does_not_reconcile_a_run_with_a_live_thread(self):
        from core.run_registry import reconcile_interrupted_runs

        run = _write_run(self.tmp, "proj", "live_run", _ORPHAN)
        touched = reconcile_interrupted_runs(runs_root=self.tmp, is_live=lambda rid: rid == "live_run")

        self.assertEqual(touched, [], "a genuinely live run must never be marked interrupted")
        self.assertFalse(self._state(run).get("interrupted"))


class LatestRunDirTests(unittest.TestCase):
    def setUp(self):
        self._prev = os.environ.get("RUNS_DIR")
        self.tmp = Path(tempfile.mkdtemp())
        os.environ["RUNS_DIR"] = str(self.tmp)
        self.project_data = {"project_name": "proj"}

    def tearDown(self):
        if self._prev is None:
            os.environ.pop("RUNS_DIR", None)
        else:
            os.environ["RUNS_DIR"] = self._prev
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_prefers_the_newest_complete_run_over_a_newer_orphan(self):
        """The false-regress guard. The newest run_id is an orphan; restore must fall to
        the older run that actually finished, not blank the pages from the orphan."""
        from ui.runtime_state import latest_run_dir

        complete = _write_run(self.tmp, "proj", "20260101T000000Z", _COMPLETE)
        _write_run(self.tmp, "proj", "20260102T000000Z", _ORPHAN)  # newer, but orphan

        self.assertEqual(latest_run_dir(self.project_data), complete)

    def test_returns_the_newest_complete_when_it_is_also_newest(self):
        from ui.runtime_state import latest_run_dir

        _write_run(self.tmp, "proj", "20260101T000000Z", _COMPLETE)
        newest = _write_run(self.tmp, "proj", "20260103T000000Z", _COMPLETE)

        self.assertEqual(latest_run_dir(self.project_data), newest)

    def test_returns_none_when_every_run_is_an_orphan(self):
        from ui.runtime_state import latest_run_dir

        _write_run(self.tmp, "proj", "20260101T000000Z", _ORPHAN)
        _write_run(self.tmp, "proj", "20260102T000000Z", _ORPHAN)

        self.assertIsNone(latest_run_dir(self.project_data), "no complete run -> None, not a blank orphan")

    def test_a_failed_but_finished_run_still_restores(self):
        """A crashed run writes done:true with its records, so its succeeded stages must
        still restore — only OOM orphans (no done) are excluded."""
        from ui.runtime_state import latest_run_dir

        failed = _write_run(
            self.tmp, "proj", "20260101T000000Z",
            {"done": True, "overall": "failure", "stages": {"QC": "success", "Analysis": "failure"}},
        )
        self.assertEqual(latest_run_dir(self.project_data), failed)

    def test_a_reconciled_orphan_is_excluded_from_restore(self):
        """Reconciliation + latest_run_dir together: after an orphan is marked
        interrupted, it stays out of the completed set that drives restore."""
        from core.run_registry import reconcile_interrupted_runs
        from ui.runtime_state import latest_run_dir

        complete = _write_run(self.tmp, "proj", "20260101T000000Z", _COMPLETE)
        _write_run(self.tmp, "proj", "20260102T000000Z", _ORPHAN)
        reconcile_interrupted_runs(runs_root=self.tmp)

        self.assertEqual(latest_run_dir(self.project_data), complete)


class RestoreStateTests(unittest.TestCase):
    """restore_run_state reads latest_run_dir; it must resurrect a complete run and
    decline when only orphans exist."""

    def setUp(self):
        self._prev = os.environ.get("RUNS_DIR")
        self.tmp = Path(tempfile.mkdtemp())
        os.environ["RUNS_DIR"] = str(self.tmp)
        self.project_data = {"project_name": "proj"}

    def tearDown(self):
        if self._prev is None:
            os.environ.pop("RUNS_DIR", None)
        else:
            os.environ["RUNS_DIR"] = self._prev
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_restores_records_from_a_complete_run(self):
        from ui.runtime_state import restore_run_state

        _write_run(
            self.tmp, "proj", "20260101T000000Z",
            {"done": True, "stages": {"QC": "success"}, "qc_record": {"status": "success"}},
        )
        restored = restore_run_state(self.project_data)
        self.assertIsNotNone(restored)
        self.assertEqual(restored["qc"]["status"], "success")

    def test_resumes_an_in_progress_run_with_config_but_no_records(self):
        """A live run (no done, NOT interrupted) must come back as config-only so a fresh
        tab reconnects and the poller keeps filling it — records stay withheld until done."""
        from ui.runtime_state import restore_run_state

        run = _write_run(self.tmp, "proj", "20260202T000000Z", _ORPHAN)
        restored = restore_run_state(self.project_data)
        self.assertIsNotNone(restored)
        self.assertEqual(restored["run_config"]["run_id"], "20260202T000000Z")
        self.assertEqual(restored["run_config"]["output_dir"], str(run))
        for stage in ("qc", "preprocess", "analysis", "databases", "phylo"):
            self.assertEqual(restored[stage], {}, f"{stage} must stay withheld mid-run")

    def test_declines_a_reconciled_interrupted_orphan(self):
        """Once reconciliation marks the orphan interrupted, restore must NOT resurrect
        it — otherwise it comes back as the active run and unlocks its blank pages."""
        from core.run_registry import reconcile_interrupted_runs
        from ui.runtime_state import restore_run_state

        _write_run(self.tmp, "proj", "20260102T000000Z", _ORPHAN)
        reconcile_interrupted_runs(runs_root=self.tmp)
        self.assertIsNone(restore_run_state(self.project_data), "an interrupted-only project restores nothing")


class StartupWiringTests(unittest.TestCase):
    def test_the_server_startup_runs_reconciliation(self):
        src = (Path(__file__).resolve().parent.parent / "app.py").read_text(encoding="utf-8")
        self.assertIn("reconcile_interrupted_runs", src, "main() must reconcile orphans at startup")


if __name__ == "__main__":
    unittest.main()
