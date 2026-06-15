"""Thread-safe pipeline progress writer."""

from __future__ import annotations

import json
import threading
from pathlib import Path

_lock = threading.Lock()


def write_progress(progress_path: Path, update: dict) -> None:
	progress_path.parent.mkdir(parents=True, exist_ok=True)
	with _lock:
		try:
			existing = json.loads(progress_path.read_text()) if progress_path.exists() else {}
		except Exception:
			existing = {}
		for key, value in update.items():
			if key == "stages" and isinstance(value, dict) and isinstance(existing.get(key), dict):
				merged = dict(existing.get(key) or {})
				merged.update(value)
				existing[key] = merged
			elif key == "log_lines" and isinstance(value, list) and isinstance(existing.get(key), list):
				existing[key] = list(existing.get(key) or []) + value
			else:
				existing[key] = value
		progress_path.write_text(json.dumps(existing), encoding="utf-8")


def append_log(progress_path: Path, entry: dict) -> None:
	progress_path.parent.mkdir(parents=True, exist_ok=True)
	with _lock:
		try:
			existing = json.loads(progress_path.read_text()) if progress_path.exists() else {}
		except Exception:
			existing = {}
		lines = existing.get("log_lines") or []
		lines.append(entry)
		existing["log_lines"] = lines
		progress_path.write_text(json.dumps(existing), encoding="utf-8")


def read_progress(progress_path: Path) -> dict:
	with _lock:
		try:
			return json.loads(progress_path.read_text()) if progress_path.exists() else {}
		except Exception:
			return {}