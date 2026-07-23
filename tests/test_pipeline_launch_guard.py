"""A pipeline may start only on an explicit click, and never twice at once.

The Setup "Run Pipeline" button used to launch the background thread itself and then
navigate to Progress, so a run — hours of remote NCBI/BOLD queries — began the instant
you left Setup, with no explicit "start now" gate. On this box (3 GB, active OOM killer)
that is dangerous, and two concurrent runs is a guaranteed OOM. The decision: the
Progress page's "Start Pipeline" button is the single, deliberate trigger, disabled
while a run is live, with the real guard in the launch function itself.

The guard keys on the LIVE THREAD, never on progress.json. This is not a nicety: at the
time of writing, runs/ already holds several runs whose progress.json reads
``stages: {"Analysis": "running"}`` with no ``done`` — runs the OOM killer took mid-flight.
A disk-based guard would read those as "still running" and refuse every future run in
those projects forever. Only the process knows what is truly executing.

Seen to fail before the fix: with Setup still calling ``_launch_pipeline_run``,
``test_setup_does_not_launch_a_run`` finds the call and goes red for the documented
reason. The OOM-safety and double-launch tests characterise the guard's contract and
would go red if it were ever rewritten against disk state.
"""

from __future__ import annotations

import json
import threading
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


class _StubRunConfig:
    """The launcher only reads run_id and output_dir; the thread body is stubbed."""

    def __init__(self, run_id: str, output_dir: Path):
        self.run_id = run_id
        self.output_dir = output_dir


class LaunchTriggerTests(unittest.TestCase):
    def test_setup_does_not_launch_a_run(self):
        """Setup configures and navigates; it must never start the thread. That is the
        auto-start this task removes — a run may begin only from the explicit button."""
        src = (ROOT / "ui" / "callbacks_setup.py").read_text(encoding="utf-8")
        self.assertNotIn(
            "_launch_pipeline_run",
            src,
            "callbacks_setup.py must not launch a run — that is the auto-start",
        )

    def test_the_launcher_is_wired_to_exactly_one_trigger(self):
        """`_launch_pipeline_run` is the sole launcher; it must be invoked from exactly
        one place — the Progress "Start Pipeline" callback — never a second trigger."""
        import subprocess

        out = subprocess.run(
            ["grep", "-rn", r"_launch_pipeline_run(", "--include=*.py", "ui/"],
            capture_output=True,
            text=True,
        )
        # Exclude the definition line; every remaining hit is a call site.
        call_sites = [
            line for line in out.stdout.splitlines() if "def _launch_pipeline_run(" not in line
        ]
        self.assertEqual(
            len(call_sites),
            1,
            "the launcher must have exactly one call site (the Start button):\n"
            + "\n".join(call_sites),
        )
        self.assertIn("callbacks_progress.py", call_sites[0])


class BackendGuardTests(unittest.TestCase):
    """The real double-launch guard lives in the launch function, not the UI: a fast
    double click or a second tab must not slip a second thread past it."""

    def setUp(self):
        import app
        import ui.callbacks_progress as cp

        self.app = app
        self.cp = cp
        self._original_thread_body = cp._run_pipeline_thread
        self._releases: list[threading.Event] = []
        app._pipeline_threads.clear()

    def tearDown(self):
        for event in self._releases:
            event.set()
        for thread in list(self.app._pipeline_threads.values()):
            thread.join(timeout=5)
        self.app._pipeline_threads.clear()
        self.cp._run_pipeline_thread = self._original_thread_body

    def _blocking_body(self):
        """A stub thread body that runs (no real tool) until released."""
        release = threading.Event()
        self._releases.append(release)
        started: list[str] = []

        def body(run_config, run_config_data, progress_path):
            started.append(run_config.run_id)
            release.wait(5)

        return body, started, release

    def test_a_second_launch_while_a_run_is_live_is_refused(self):
        body, started, _ = self._blocking_body()
        self.cp._run_pipeline_thread = body
        with __import__("tempfile").TemporaryDirectory() as d:
            rc1 = _StubRunConfig("runA", Path(d) / "a")
            rc2 = _StubRunConfig("runB", Path(d) / "b")

            self.assertTrue(self.cp._launch_pipeline_run(rc1, {}), "first run must start")
            # A different run_id, but a run is already live: refused. On this box two
            # pipelines at once is OOM, so the guard is global, not per-run_id.
            self.assertFalse(self.cp._launch_pipeline_run(rc2, {}), "second run must be refused")
            self.assertEqual(started, ["runA"], "the refused run must never have started a thread")

    def test_a_stale_running_progress_json_without_a_thread_does_not_block_launch(self):
        """OOM-safety. A run killed mid-flight leaves stages=running / no done on disk.
        The guard must ignore that and launch, because no thread is alive."""
        done = threading.Event()

        def body(run_config, run_config_data, progress_path):
            done.set()

        self.cp._run_pipeline_thread = body
        with __import__("tempfile").TemporaryDirectory() as d:
            out_dir = Path(d) / "orphan"
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "progress.json").write_text(
                json.dumps({"run_id": "orphan", "stages": {"Analysis": "running"}}),
                encoding="utf-8",
            )
            # No thread registered — this is the OOM orphan case.
            self.assertFalse(self.cp._any_pipeline_running())
            rc = _StubRunConfig("orphan", out_dir)
            self.assertTrue(
                self.cp._launch_pipeline_run(rc, {}),
                "a stale disk 'running' with no live thread must not block a fresh launch",
            )
            self.assertTrue(done.wait(5), "the launched thread should have run")

    def test_any_pipeline_running_tracks_the_live_thread_not_disk(self):
        """The Start button reads this, so it must reflect the process, not the file."""
        self.assertFalse(self.cp._any_pipeline_running(), "nothing running yet")

        body, _started, release = self._blocking_body()
        self.cp._run_pipeline_thread = body
        with __import__("tempfile").TemporaryDirectory() as d:
            rc = _StubRunConfig("live", Path(d) / "live")
            self.assertTrue(self.cp._launch_pipeline_run(rc, {}))
            self.assertTrue(self.cp._any_pipeline_running(), "a live thread must read as running")
            release.set()
            for thread in list(self.app._pipeline_threads.values()):
                thread.join(timeout=5)
        self.assertFalse(self.cp._any_pipeline_running(), "after the thread ends, not running")


class ButtonDisabledWiringTests(unittest.TestCase):
    """The button's disabled state must be wired to the live-thread signal in both the
    click and the poll callbacks — the click for immediacy, the poll to re-enable when
    the run ends (or is OOM-killed)."""

    def test_both_callbacks_drive_the_start_button_disabled_from_the_live_signal(self):
        src = (ROOT / "ui" / "callbacks_progress.py").read_text(encoding="utf-8")
        self.assertEqual(
            src.count('Output("btn-start-progress", "disabled"'),
            2,
            "both _start_pipeline and _poll_progress must output btn-start-progress.disabled",
        )
        # And that output is fed by the live-thread helper, never by a disk 'running'.
        self.assertIn("_any_pipeline_running()", src)


if __name__ == "__main__":
    unittest.main()
