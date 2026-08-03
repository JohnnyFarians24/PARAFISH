"""Task 2: a stage record belongs to the run that produced it, and only to that run.

The stage stores are per-session and the run-start callback writes only
``store-run-config``. Nothing cleared the five ``store-*-outputs``, so starting a new
run left the previous run's records in place with ``status: "success"`` — every stage
page reads the store before the disk, so the tabs stayed unlocked on the OLD run's
tables while a new run was executing.

The measured mechanism was NOT ``latest_run_dir``: ``_resolve_output_dir`` already
prefers the active run's ``output_dir`` and correctly returned the new run. It was the
stale store, consulted first.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from ui.pipeline_state import RUN_ID_KEY, record_belongs_to, stamp_record


def _output_key(client, output_id: str, input_id: str) -> tuple[str, list]:
    """The dependency key Dash registered for this callback.

    ``allow_duplicate`` outputs carry an ``@hash`` suffix, so the key cannot be
    rebuilt by hand from the component ids.
    """
    for dep in client.get("/_dash-dependencies").get_json():
        if output_id in dep["output"] and any(i["id"] == input_id for i in dep["inputs"]):
            outputs = [
                {"id": o.split(".")[0], "property": o.split(".")[1].split("@")[0]}
                for o in dep["output"].strip(".").split("...")
            ]
            return dep["output"], outputs, dep["inputs"]
    raise AssertionError(f"no callback outputs {output_id} from {input_id}")


def _write_run(root: Path, run_id: str, *, done: bool, records: bool = True) -> Path:
    run_dir = root / run_id
    (run_dir / "databases").mkdir(parents=True)
    state = {"run_id": run_id, "done": done, "stages": {"databases": "success" if done else "running"}}
    if records:
        for key in ("qc_record", "preprocess_record", "analysis_record", "db_record", "phylo_record"):
            state[key] = {"status": "success", "outputs": {"s1": {}}}
    (run_dir / "progress.json").write_text(json.dumps(state), encoding="utf-8")
    return run_dir


class RecordIdentityTests(unittest.TestCase):
    def test_an_unstamped_record_never_matches(self):
        """Records already sitting in open sessions carry no stamp. Treating them as
        compatible would keep the bug alive for exactly the users who have the app
        open — the population this fix is for."""
        self.assertFalse(record_belongs_to({"status": "success"}, "20260724T000000Z"))

    def test_a_record_matches_only_its_own_run(self):
        record = stamp_record({"status": "success"}, "20260724T000000Z")

        self.assertTrue(record_belongs_to(record, "20260724T000000Z"))
        self.assertFalse(record_belongs_to(record, "20260101T000000Z"))

    def test_stamping_an_empty_record_stays_empty(self):
        """An absent record must not become a truthy dict carrying only a stamp —
        `if any(bool(store) ...)` in the restore gate keys on exactly that."""
        self.assertEqual(stamp_record({}, "r1"), {})
        self.assertEqual(stamp_record(None, "r1"), {})

    def test_stamp_does_not_disturb_the_record(self):
        record = stamp_record({"status": "success", "outputs": {"s1": {}}}, "r1")

        self.assertEqual(record["status"], "success")
        self.assertEqual(record["outputs"], {"s1": {}})
        self.assertEqual(record[RUN_ID_KEY], "r1")


class RestoreStampConsistencyTests(unittest.TestCase):
    """The condition this task carries the most risk on.

    ``restore_run_state`` feeds the four Eixo B pages on a fresh session. If it stamped
    with a run_id that disagreed with the run_config it returns alongside, the identity
    check would reject every restored record and blank all four pages at once — the
    restore bug, reintroduced by its own guard.
    """

    def _project(self, root: Path) -> dict:
        return {"active_id": "p1", "projects": {"p1": {"project_name": root.name}}}

    def test_restored_records_match_the_run_config_restored_with_them(self):
        from ui.runtime_state import restore_run_state

        with tempfile.TemporaryDirectory() as d:
            runs_root = Path(d) / "runs"
            project = runs_root / "proj"
            project.mkdir(parents=True)
            _write_run(project, "20260724T000000Z", done=True)
            with patch("ui.runtime_state.runs_root", return_value=runs_root):
                restored = restore_run_state(self._project(project))

        self.assertIsNotNone(restored)
        run_id = restored["run_config"]["run_id"]
        for stage in ("qc", "preprocess", "analysis", "databases", "phylo"):
            with self.subTest(stage=stage):
                self.assertTrue(
                    record_belongs_to(restored[stage], run_id),
                    f"{stage} restored but would be rejected as another run's record",
                )

    def test_a_live_run_restores_the_config_without_unlocking_its_stages(self):
        """Pre-existing behaviour that must survive: a run still going hands back its
        run_config (so the poller resumes) but no stage records."""
        from ui.runtime_state import restore_run_state

        with tempfile.TemporaryDirectory() as d:
            runs_root = Path(d) / "runs"
            project = runs_root / "proj"
            project.mkdir(parents=True)
            _write_run(project, "20260724T000000Z", done=False)
            with patch("ui.runtime_state.runs_root", return_value=runs_root):
                restored = restore_run_state(self._project(project))

        self.assertEqual(restored["run_config"]["run_id"], "20260724T000000Z")
        for stage in ("qc", "preprocess", "analysis", "databases", "phylo"):
            self.assertEqual(restored[stage], {})


class StaleRecordDroppedTests(unittest.TestCase):
    STORES = ["qc", "preprocess", "analysis", "databases", "phylo"]

    def _fire(self, run_config, store_value):
        from app import create_app

        client = create_app().server.test_client()
        key, outputs, _ = _output_key(client, "store-qc-outputs", "store-run-config")
        body = {
            "output": key,
            "outputs": outputs,
            "inputs": [{"id": "store-run-config", "property": "data", "value": run_config}],
            "state": [
                {"id": f"store-{s}-outputs", "property": "data", "value": store_value}
                for s in self.STORES
            ],
            "changedPropIds": ["store-run-config.data"],
        }
        return client.post("/_dash-update-component", json=body)

    def test_the_previous_runs_records_do_not_survive_into_a_new_run(self):
        """The reported symptom, end to end through the real callback."""
        old = stamp_record({"status": "success", "outputs": {"s1": {}}}, "20260101T000000Z")
        response = self._fire({"run_id": "20260724T000000Z", "output_dir": "/tmp/x"}, old)

        payload = response.get_json()["response"]
        for stage in self.STORES:
            with self.subTest(stage=stage):
                self.assertEqual(payload[f"store-{stage}-outputs"]["data"], {})

    def test_an_unstamped_record_is_dropped_too(self):
        """Sessions open when this shipped hold unstamped records. Strict by design."""
        response = self._fire({"run_id": "r-new", "output_dir": "/tmp/x"}, {"status": "success"})

        payload = response.get_json()["response"]
        self.assertEqual(payload["store-qc-outputs"]["data"], {})

    def test_the_current_runs_records_are_left_alone(self):
        """Not a clear-on-start: a matching record survives re-renders. Every store
        unchanged means PreventUpdate, which Dash answers with an empty 204."""
        mine = stamp_record({"status": "success"}, "20260724T000000Z")
        response = self._fire({"run_id": "20260724T000000Z", "output_dir": "/tmp/x"}, mine)

        self.assertFalse(response.get_json())


class WaitingBannerTests(unittest.TestCase):
    STAGES = ["qc", "preprocess", "analysis", "databases", "phylo"]

    def _banners(self, run_config):
        from app import create_app

        client = create_app().server.test_client()
        key, outputs, inputs = _output_key(client, "databases-run-banner", "store-run-config")
        values = {"store-run-config": run_config}
        body = {
            "output": key,
            "outputs": outputs,
            "inputs": [
                {"id": i["id"], "property": i["property"], "value": values.get(i["id"], {})}
                for i in inputs
            ],
            "state": [],
            "changedPropIds": ["store-run-config.data"],
        }
        return client.post("/_dash-update-component", json=body).get_json()["response"]

    def test_a_stage_pending_in_the_active_run_says_so(self):
        with tempfile.TemporaryDirectory() as d:
            run_dir = _write_run(Path(d), "20260724T000000Z", done=False, records=False)
            payload = self._banners({"run_id": "20260724T000000Z", "output_dir": str(run_dir)})

        banner = json.dumps(payload["databases-run-banner"]["children"], ensure_ascii=False)
        self.assertIn("has not run in this run yet", banner)

    def test_no_banner_when_the_run_has_finished(self):
        """A finished run has nothing pending; the page's own empty state stands."""
        with tempfile.TemporaryDirectory() as d:
            run_dir = _write_run(Path(d), "20260724T000000Z", done=True)
            payload = self._banners({"run_id": "20260724T000000Z", "output_dir": str(run_dir)})

        self.assertEqual(payload["databases-run-banner"]["children"], "")

    def test_no_banner_without_an_active_run(self):
        self.assertEqual(self._banners({})["databases-run-banner"]["children"], "")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
