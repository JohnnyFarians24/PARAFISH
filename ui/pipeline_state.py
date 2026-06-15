"""Shared helpers for storing and reading pipeline stage state in Dash stores."""

from __future__ import annotations

import json
from copy import deepcopy
from datetime import datetime, timezone
from typing import Any, Optional


def utc_now_iso() -> str:
	return datetime.now().astimezone().replace(microsecond=0).isoformat()


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
