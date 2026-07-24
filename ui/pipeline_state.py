"""Shared helpers for storing and reading pipeline stage state in Dash stores."""

from __future__ import annotations

import json
from copy import deepcopy
from datetime import datetime, timezone
from typing import Any, Optional


def utc_now_iso() -> str:
	return datetime.now().astimezone().replace(microsecond=0).isoformat()


# The run a stage record was produced by. Dunder-ish so it cannot collide with a key a
# StepOutput record legitimately carries (status, outputs, metrics_files, ...).
RUN_ID_KEY = "__run_id__"


def stamp_record(record: Optional[dict], run_id: Optional[str]) -> dict:
	"""Tag ``record`` with the run that produced it.

	Every write to a stage store goes through here. The stage stores are per-session
	and nothing cleared them when a new run started, so the previous run's records
	survived into the new run with ``status: "success"`` intact — every stage page
	read them in preference to the disk, and showed the old run's results while a new
	one was executing. A record that cannot say which run it came from cannot be
	checked against the run on screen.
	"""
	stamped = dict(record or {})
	if not stamped:
		return stamped
	stamped[RUN_ID_KEY] = str(run_id or "")
	return stamped


def record_belongs_to(record: Optional[dict], run_id: Optional[str]) -> bool:
	"""Whether ``record`` was produced by ``run_id``.

	An UNSTAMPED record is never a match. Records written before this existed are
	already sitting in open sessions, and treating them as compatible would keep the
	bug alive for exactly the users who have the app open — which is the population
	the fix is for. The cost of the strict reading is one refresh; the cost of the
	lenient one is trusting data that cannot identify itself, which is the failure
	this whole mechanism exists to close.
	"""
	stamped = str((record or {}).get(RUN_ID_KEY) or "")
	return bool(stamped) and stamped == str(run_id or "")


def get_stage_record(step_outputs: Optional[dict], stage: str) -> dict:
	return dict((step_outputs or {}).get(stage) or {})


def merge_stage_record(step_outputs: Optional[dict], stage: str, record: dict) -> dict:
	merged = deepcopy(step_outputs or {})
	merged[stage] = record
	return merged


def _normalize_store(data: Optional[dict], _depth: int = 0) -> dict:
	"""Normalize a store payload to a dict of dict values."""
	if _depth > 5:
		return data if isinstance(data, dict) else {}
	if not data:
		return {}
	if isinstance(data, str):
		try:
			data = json.loads(data)
		except Exception:
			return {}
	if not isinstance(data, dict):
		return {}
	result = {}
	for key, value in data.items():
		if isinstance(value, dict):
			result[key] = _normalize_store(value, _depth + 1)
		elif isinstance(value, str):
			try:
				parsed = json.loads(value)
			except Exception:
				result[key] = value
			else:
				if isinstance(parsed, dict):
					result[key] = _normalize_store(parsed, _depth + 1)
				elif isinstance(parsed, list):
					result[key] = parsed
				else:
					result[key] = value
		else:
			result[key] = value
	return result


def get_stage_outputs(step_outputs: Optional[dict], stage: str) -> dict:
	record = get_stage_record(step_outputs, stage)
	outputs = record.get("outputs") or {}
	return outputs if isinstance(outputs, dict) else {}


def get_stage_error(step_outputs: Optional[dict], stage: str) -> str:
	record = get_stage_record(step_outputs, stage)
	return str(record.get("error") or record.get("error_message") or "")


def stage_status(record: Optional[dict]) -> str:
	return str((record or {}).get("status") or "pending")


def stage_ready(record: Optional[dict]) -> bool:
	return stage_status(record) == "success"


def stage_failed(record: Optional[dict]) -> bool:
	return stage_status(record) in {"failure", "blocked"}
