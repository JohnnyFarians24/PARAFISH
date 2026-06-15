"""Path management utilities for projects, runs and step outputs."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from datetime import datetime


UPLOAD_DIRECTORY = os.path.join(os.getcwd(), "uploads")
os.makedirs(UPLOAD_DIRECTORY, exist_ok=True)

PROJECTS_ROOT = os.path.join(os.getcwd(), "projects")
os.makedirs(PROJECTS_ROOT, exist_ok=True)


@dataclass
class ProjectPaths:
	"""Resolved path references used by pipeline stages."""

	root: str
	data: str
	outputs: str
	logs: str


def _safe_name(text):
	value = str(text or "").strip()
	if not value:
		return "Project"
	cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
	cleaned = cleaned.strip("._")
	return cleaned or "Project"


def _project_dir(project_name):
	return os.path.join(PROJECTS_ROOT, _safe_name(project_name))


def _project_exports_dir(project_dir):
	if not project_dir:
		return ""
	exports_dir = os.path.join(project_dir, "exports")
	os.makedirs(exports_dir, exist_ok=True)
	return exports_dir


def _ensure_project(project_name):
	safe_project = _safe_name(project_name)
	pdir = _project_dir(safe_project)
	os.makedirs(pdir, exist_ok=True)
	os.makedirs(os.path.join(pdir, "runs"), exist_ok=True)
	os.makedirs(os.path.join(pdir, "uploads"), exist_ok=True)
	os.makedirs(_project_exports_dir(pdir), exist_ok=True)
	meta_path = os.path.join(pdir, "project.json")
	if not os.path.exists(meta_path):
		from .utils import _write_json

		_write_json(
			meta_path,
			{
				"project_name": safe_project,
				"created_at": datetime.now().isoformat(timespec="seconds"),
			},
		)
	return pdir


def _create_run(project_dir):
	runs_dir = os.path.join(project_dir, "runs")
	os.makedirs(runs_dir, exist_ok=True)
	stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
	run_dir = os.path.join(runs_dir, stamp)
	os.makedirs(run_dir, exist_ok=True)
	return run_dir


def _get_active_root(project_state):
	state = project_state or {}
	run_dir = state.get("run_dir")
	if run_dir:
		return run_dir
	project_dir = state.get("project_dir")
	if project_dir:
		return project_dir
	return UPLOAD_DIRECTORY


def _step_base(project_state, step_name):
	base = _get_active_root(project_state)
	if base == UPLOAD_DIRECTORY:
		return base
	step_dir = os.path.join(base, "steps", step_name)
	os.makedirs(step_dir, exist_ok=True)
	return step_dir


def _step_outputs_base(project_state, step_name):
	step_dir = _step_base(project_state, step_name)
	outputs_dir = os.path.join(step_dir, "outputs")
	os.makedirs(outputs_dir, exist_ok=True)
	return outputs_dir


def _step_output_path(project_state, step_name, filename, sample_id=None):
	outputs_dir = _step_outputs_base(project_state, step_name)
	if sample_id:
		outputs_dir = os.path.join(outputs_dir, _safe_name(sample_id))
		os.makedirs(outputs_dir, exist_ok=True)
	return os.path.join(outputs_dir, filename) if filename else outputs_dir


def _resolve_output_path(project_state, step_name, output_type, filename=None, sample_id=None):
	base = _get_active_root(project_state)
	if base == UPLOAD_DIRECTORY:
		base = project_state.get("run_dir") or project_state.get("project_dir") or UPLOAD_DIRECTORY

	steps_root = os.path.join(base, "steps")
	results_root = os.path.join(base, "results")
	os.makedirs(steps_root, exist_ok=True)
	os.makedirs(results_root, exist_ok=True)

	if output_type == "intermediate":
		return _step_output_path(project_state, step_name, filename, sample_id=sample_id)

	if output_type == "final":
		final_dir = os.path.join(results_root, step_name)
		if sample_id:
			final_dir = os.path.join(final_dir, _safe_name(sample_id))
		os.makedirs(final_dir, exist_ok=True)
		return os.path.join(final_dir, filename) if filename else final_dir

	raise ValueError(f"Unknown output_type: {output_type}")


def _project_session_path(project_dir):
	return os.path.join(project_dir, "session_state.json")


def _run_session_path(run_dir):
	return os.path.join(run_dir, "session.json")
