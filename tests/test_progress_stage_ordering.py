"""Stage badges must not lag a step behind the pipeline.

The Progress page showed "Analysis: running" while the log had already written
"Starting Databases stage": every completed stage stayed "running" on disk until the
NEXT stage's write finally flipped it. Two independent causes, both in
``_run_pipeline_thread``:

  * the completing stage set ``stages[X] = record["status"]`` AFTER the write that
    persisted both the record and the badge, so a snapshot could hold a *success
    record* next to a *running badge* — and that terminal status only reached disk on
    the next stage's write; and
  * ``_publish()`` was a no-op, so a stage's "running" never reached disk when it
    started — only when it ended.

This is the first UI task whose file — ``progress.json`` — the six-layer golden does
NOT freeze, so ``golden.py check all`` stays green even if the state writer breaks.
The check therefore has to carry the whole weight itself: it drives the REAL
``_run_pipeline_thread`` (QC+Preprocess+Analysis enabled, zero samples so no external
tool runs), snapshots progress.json after every write, and proves three things that
the golden cannot:

  1. ORDERING — no on-disk snapshot ever holds a terminal record beside a "running"
     badge, and each stage is already terminal when the next is announced;
  2. COMPLETENESS — the finished progress.json still carries run_id, run_config, all
     five stage records, log_lines, stages, done and overall; and
  3. RESTORE — ``restore_run_state`` (the disk-reconstruction that reads progress.json
     in a fresh tab) still rebuilds every stage record from the produced file.

Seen to fail before the fix: with the pre-fix ordering, the QC snapshot carries
``qc_record.status == "success"`` while ``stages["QC"] == "running"`` — assertion (1)
goes red for the documented reason. (Demonstrated by reverting callbacks_progress.py
and re-running: the ORDERING test fails, the other two still pass.)
"""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

TERMINAL = {"success", "failure", "skipped"}
# Badge key in the "stages" map -> the record key that carries the same stage's status.
STAGE_RECORD = {
    "QC": "qc_record",
    "Preprocess": "preprocess_record",
    "Analysis": "analysis_record",
    "Databases": "db_record",
    "Phylo": "phylo_record",
}


def _load_run_pipeline():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    spec = importlib.util.spec_from_file_location("run_pipeline", ROOT / "scripts" / "run_pipeline.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _Driven:
    """The outcome of one real, offline pipeline run."""

    def __init__(self, snapshots, final, project_name, run_dir):
        self.snapshots = snapshots
        self.final = final
        self.project_name = project_name
        self.run_dir = run_dir


def _drive(tmp: Path, enabled: list) -> _Driven:
    """Run the real _run_pipeline_thread with no samples, capturing every disk write."""
    run_pipeline = _load_run_pipeline()
    import ui.callbacks_progress as cp
    import ui.runtime_state as rs
    from core.progress_writer import append_log as real_append
    from core.progress_writer import read_progress
    from core.progress_writer import write_progress as real_write
    from ui.runtime_state import run_config_from_store, slugify

    project_name = "ordering-probe"
    run_id = "20260101T000000Z"
    runs_root = tmp / "runs"
    out_dir = runs_root / slugify(project_name) / run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    progress_path = out_dir / "progress.json"

    project_data, samples_data, run_config_data = run_pipeline.build_stores(
        run_pipeline.FIXTURE_INPUTS, out_dir, enabled, ["silva"], run_id, {"threads": 1}
    )
    run_config = run_config_from_store(project_data, samples_data, run_config_data)
    run_config.run_id = run_id
    run_config.output_dir = out_dir
    # No samples: the per-sample loops run no external tool, but every stage still
    # executes its running -> terminal -> write path, which is what we are probing.
    run_config.samples = []

    snapshots: list[dict] = []

    def _snap() -> None:
        snapshots.append(copy.deepcopy(read_progress(progress_path)))

    def _wrapped_write(path, update):
        real_write(path, update)
        _snap()

    def _wrapped_append(path, entry):
        real_append(path, entry)
        _snap()

    original_runs_root = rs.runs_root
    cp.write_progress = _wrapped_write
    cp.append_log = _wrapped_append
    rs.runs_root = lambda: runs_root
    try:
        cp._run_pipeline_thread(run_config, run_config_data, progress_path)
    finally:
        cp.write_progress = real_write
        cp.append_log = real_append
        rs.runs_root = original_runs_root

    final = json.loads(progress_path.read_text(encoding="utf-8"))
    return _Driven(snapshots, final, project_name, out_dir)


class StageOrderingTests(unittest.TestCase):
    ENABLED = ["qc", "preprocess", "analysis"]

    def test_no_snapshot_shows_a_terminal_record_beside_a_running_badge(self):
        """The core regression. A stage whose record already reads success/failure must
        never still read "running" in the badge map of the same on-disk snapshot."""
        with tempfile.TemporaryDirectory() as d:
            driven = _drive(Path(d), self.ENABLED)

        self.assertTrue(driven.snapshots, "the run produced no progress.json writes")
        violations = []
        for i, snap in enumerate(driven.snapshots):
            stages = snap.get("stages") or {}
            for badge, rec_key in STAGE_RECORD.items():
                status = (snap.get(rec_key) or {}).get("status")
                if status in TERMINAL and stages.get(badge) == "running":
                    violations.append(f"snapshot {i}: {rec_key}={status} but stages[{badge}]=running")
        self.assertEqual(violations, [], "badge lagged its record:\n" + "\n".join(violations))

    def test_each_stage_is_terminal_when_the_next_is_announced(self):
        """At the instant "Starting <N+1> stage" is logged, stage <N>'s badge must
        already be terminal — the exact symptom from the video."""
        with tempfile.TemporaryDirectory() as d:
            driven = _drive(Path(d), self.ENABLED)

        prev_of = {"Preprocess": "QC", "Analysis": "Preprocess"}
        checked = set()
        for snap in driven.snapshots:
            lines = snap.get("log_lines") or []
            if not lines:
                continue
            msg = lines[-1].get("msg", "")
            for stage, prev in prev_of.items():
                if msg == f"Starting {stage} stage":
                    badge = (snap.get("stages") or {}).get(prev)
                    self.assertIn(
                        badge,
                        TERMINAL,
                        f"{prev} badge was {badge!r} when {stage} was announced",
                    )
                    checked.add(stage)
        # Prove the boundaries actually occurred, so a run that silently skipped a
        # stage cannot pass this test vacuously.
        self.assertEqual(checked, set(prev_of), f"missing stage boundaries: {set(prev_of) - checked}")

    def test_a_stage_shows_running_on_disk_the_moment_it_starts(self):
        """`_publish()` used to be a no-op, so "running" only reached disk when the
        stage ended. Every enabled stage must appear "running" in some snapshot."""
        with tempfile.TemporaryDirectory() as d:
            driven = _drive(Path(d), self.ENABLED)

        for badge in ("QC", "Preprocess", "Analysis"):
            ran = any((snap.get("stages") or {}).get(badge) == "running" for snap in driven.snapshots)
            self.assertTrue(ran, f"{badge} never appeared as 'running' on disk")


class ProgressFileCompletenessTests(unittest.TestCase):
    """The golden does not freeze progress.json, so the reorder must be proven not to
    have dropped anything the file has always carried."""

    def test_finished_progress_json_still_carries_every_key(self):
        with tempfile.TemporaryDirectory() as d:
            driven = _drive(Path(d), ["qc", "preprocess", "analysis"])
        final = driven.final

        self.assertEqual(final.get("run_id"), "20260101T000000Z")
        self.assertTrue(final.get("done"))
        self.assertIn(final.get("overall"), ("success", "failure"))
        self.assertIsInstance(final.get("run_config"), dict)
        self.assertIn("db_sources", final["run_config"])
        self.assertIsInstance(final.get("stages"), dict)
        self.assertTrue(final.get("log_lines"), "log_lines must survive the run")
        for rec_key in STAGE_RECORD.values():
            record = final.get(rec_key)
            self.assertIsInstance(record, dict, f"{rec_key} missing")
            self.assertIn("status", record, f"{rec_key} has no status")


class ProgressRestoreTests(unittest.TestCase):
    """restore_run_state reads progress.json to rebuild a fresh tab. The reorder must
    not disturb the records it hands back."""

    def test_restore_rebuilds_every_stage_record_from_the_finished_run(self):
        with tempfile.TemporaryDirectory() as d:
            driven = _drive(Path(d), ["qc", "preprocess", "analysis"])
            import ui.runtime_state as rs

            original = rs.runs_root
            rs.runs_root = lambda: driven.run_dir.parent.parent
            try:
                restored = rs.restore_run_state({"project_name": driven.project_name})
            finally:
                rs.runs_root = original

        self.assertIsNotNone(restored, "a finished run must be restorable")
        self.assertEqual(restored["run_config"]["run_id"], "20260101T000000Z")
        # done is True, so every stage record is handed back with its status.
        for stage in ("qc", "preprocess", "analysis"):
            self.assertIn("status", restored[stage], f"{stage} record not restored")


if __name__ == "__main__":
    unittest.main()
