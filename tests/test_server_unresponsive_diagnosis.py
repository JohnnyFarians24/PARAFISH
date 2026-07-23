"""Task 4 — "Callback failed: the server did not respond".

DIAGNOSIS (measured, not assumed):

  (a) OOM / memory pressure — CONFIRMED and DOMINANT, environmental. dmesg on this box
      carries `Out of memory: Killed process <pid> (MainThread)`; free memory sat at
      ~38 MB of 2840 with swap exhausted. When the OS SIGKILLs the Flask process (or it
      thrashes on swap) the server stops answering and every in-flight callback reports
      "the server did not respond". This is the machine, not the code.

  (b) Expensive polling tick — MEASURED CHEAP, not the cause. A full real-run
      progress.json is 33 KB with 28 log lines and parses in ~0.12 ms per tick. Re-reading
      it every 2 s costs nothing and frees immediately; there is no large-file re-parse to
      optimise. Inventing an "incremental read" here would be code pretending to fix a
      machine problem.

  (c) Heavy render during the run — ALREADY MITIGATED. The progress poll withholds all
      five stage-record store outputs (`no_update`) until the whole run is `done`, so the
      heavy per-stage tables and plots never rebuild mid-run.

CONCLUSION: the dominant cause is environmental (RAM). Per the task, it is reported as a
system requirement and NOT worked around in code. The honest user-facing surface already
exists and is correctly placed: the Setup preflight panel (free RAM, recommended sample
cap, OOM warnings) and the connection banner that collapses the repeated toasts into one
line naming OOM as the cause and pointing at that panel.

This file adds no fix. It pins the two memory-relevant properties that keep the situation
from getting worse, so a future change cannot silently undo them:
  * the poll must keep withholding heavy renders until `done`; and
  * the OOM warning surfaces must stay mounted and honest.
"""

from __future__ import annotations

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


class RenderWithholdingDuringRunTests(unittest.TestCase):
    def test_the_poll_withholds_every_stage_store_until_the_run_is_done(self):
        """During the run the poll returns no_update for the five stage stores, so the
        heavy per-stage tables/plots do not rebuild every tick. If a future change
        publishes records incrementally, mid-run memory pressure returns — and OOM with
        it. Keep the withhold."""
        src = (ROOT / "ui" / "callbacks_progress.py").read_text(encoding="utf-8")
        marker = "if not done:"
        self.assertIn(marker, src, "the poll must have a not-done early return")
        window = src[src.index(marker) : src.index(marker) + 500]
        self.assertIn(
            "no_update, no_update, no_update, no_update, no_update",
            window,
            "the not-done poll path must withhold all five stage stores",
        )


class HonestOOMMessagingTests(unittest.TestCase):
    """The task's actionable half for an environmental cause: confirm the warning is
    visible in the right place, and protect it from quietly disappearing."""

    def test_the_connection_banner_is_mounted_once_in_the_app_layout(self):
        layout_src = (ROOT / "ui" / "layout.py").read_text(encoding="utf-8")
        self.assertEqual(
            layout_src.count("connection_banner()"),
            1,
            "the connection banner must be mounted exactly once in the app layout",
        )

    def test_the_banner_names_the_memory_cause_not_the_callback(self):
        """A memory-starved box's real cause is almost never the callback in flight, so
        the banner must say OOM/memory — otherwise users chase the wrong thing."""
        banner_src = (ROOT / "ui" / "components" / "connection_banner.py").read_text(encoding="utf-8")
        self.assertTrue(
            "OOM" in banner_src or "memoria" in banner_src or "memória" in banner_src,
            "the banner must name the memory (OOM) cause",
        )

    def test_the_setup_preflight_reports_free_memory(self):
        """The banner sends users to the Setup resource panel; that panel must actually
        read available memory and surface it, not just claim to."""
        setup_src = (ROOT / "ui" / "callbacks_setup.py").read_text(encoding="utf-8")
        self.assertIn("available_memory_bytes", setup_src)
        self.assertIn("RAM", setup_src)

    def test_preflight_ties_the_symptom_to_the_cause_for_maintainers(self):
        """core.preflight documents the OOM -> 'the server did not respond' chain, so the
        next maintainer does not re-diagnose it from scratch."""
        preflight_src = (ROOT / "core" / "preflight.py").read_text(encoding="utf-8")
        self.assertIn("did not respond", preflight_src)
        self.assertIn("OOM", preflight_src)


if __name__ == "__main__":
    unittest.main()
