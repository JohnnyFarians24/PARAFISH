"""Shared helpers for serializing UI state into runner-friendly runtime objects."""

from __future__ import annotations

import base64
import re
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from core.config import get_runs_dir
from core.config import get_db_paths, get_feature_flags
from core.models import AnalysisBranch, AnalysisMode, Marker, PhyloMode, RunConfig, Sample, SamplePaths, StepOutput


def slugify(value: str) -> str:
	value = re.sub(r"[^A-Za-z0-9]+", "_", value or "").strip("_")
	return value.lower() or "project"


def timestamp_slug() -> str:
	return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

def runs_root() -> Path:
	preferred = get_runs_dir()
	try:
		preferred.mkdir(parents=True, exist_ok=True)
		return preferred
	except Exception:
		fallback = Path.cwd() / "runs"
		fallback.mkdir(parents=True, exist_ok=True)
		return fallback


def project_run_dir(project_data: Optional[dict], run_id: Optional[str] = None) -> Path:
	project_data = resolve_project_data(project_data)
	project_name = project_data.get("project_name") or "project"
	effective_run_id = run_id or timestamp_slug()
	return runs_root() / slugify(project_name) / effective_run_id


def latest_run_dir(project_data: Optional[dict]) -> Optional[Path]:
	"""Newest run directory of the active project, or ``None``.

	The stage stores are per-session, so a restarted app (or a fresh tab) loses
	every result even though the run is still on disk. This is how the UI finds
	that run again.
	"""
	project = resolve_project_data(project_data)
	project_name = project.get("project_name")
	if not project_name:
		return None
	project_dir = runs_root() / slugify(project_name)
	if not project_dir.is_dir():
		return None
	runs = [path for path in project_dir.iterdir() if path.is_dir() and (path / "progress.json").exists()]
	if not runs:
		return None
	return max(runs, key=lambda path: path.stat().st_mtime)


def restore_run_state(project_data: Optional[dict]) -> Optional[dict]:
	"""Rebuild the run-config and stage records of the project's newest run.

	Returns ``None`` when there is nothing to restore. Keys mirror the session
	stores so a callback can hand them straight back to Dash.
	"""
	import json

	run_dir = latest_run_dir(project_data)
	if run_dir is None:
		return None
	try:
		state = json.loads((run_dir / "progress.json").read_text(encoding="utf-8"))
	except (OSError, ValueError):
		return None
	if not state.get("run_id"):
		return None
	# Only a FINISHED run may hand its stage records back to the stores. While a run
	# is still going, progress.json already holds the records of the stages that have
	# completed so far — restoring those here (in a fresh session/tab, where the stores
	# are empty) would unlock QC/Preprocess/Analysis one at a time, mid-run, exactly
	# what the progress poller deliberately withholds until `done`. Restore the
	# run-config regardless so the poller resumes and fills every stage together at the
	# end; leave the stage stores empty until then so the tabs stay locked.
	done = bool(state.get("done"))
	return {
		"run_config": {"run_id": state["run_id"], "output_dir": str(run_dir)},
		"qc": (state.get("qc_record") or {}) if done else {},
		"preprocess": (state.get("preprocess_record") or {}) if done else {},
		"analysis": (state.get("analysis_record") or {}) if done else {},
		"databases": (state.get("db_record") or {}) if done else {},
		"phylo": (state.get("phylo_record") or {}) if done else {},
	}


def resolve_project_data(project_store: Optional[dict]) -> dict:
	project_store = project_store or {}
	projects = project_store.get("projects")
	active_id = project_store.get("active_id")
	if isinstance(projects, dict) and active_id and active_id in projects:
		active_project = projects.get(active_id) or {}
		if isinstance(active_project, dict):
			return active_project
	if project_store.get("project_name"):
		return project_store
	return {}


def decode_upload_content(content: str) -> bytes:
	_, encoded = content.split(",", 1)
	return base64.b64decode(encoded)


def normalize_upload_filename(filename: str) -> str:
	return Path(filename).name


def enum_from_value(enum_cls, value, default):
	for item in enum_cls:
		if item.value == value or item.name == value:
			return item
	return default


def sample_record_to_sample(record: dict) -> Sample:
	mode = enum_from_value(AnalysisMode, record.get("mode"), AnalysisMode.SE)
	marker = enum_from_value(Marker, record.get("marker"), Marker.OTHER)
	paths_dict = record.get("paths") or {}
	r1_value = record.get("r1_path") or record.get("r1") or paths_dict.get("r1")
	r2_value = record.get("r2_path") or record.get("r2") or paths_dict.get("r2")
	if not r1_value:
		raise KeyError("r1_path")
	r1 = Path(r1_value)
	r2 = Path(r2_value) if r2_value else None
	group = record.get("group")
	return Sample(
		internal_id=str(record.get("internal_id") or record.get("internal_sample_id") or slugify(record.get("display_label", "sample"))),
		display_label=str(record.get("display_label") or record.get("display_name") or record.get("internal_id") or "sample"),
		paths=SamplePaths(r1=r1, r2=r2),
		mode=mode,
		marker=marker,
		selected=bool(record.get("selected", True)),
		is_control=bool(record.get("is_control", False)),
		group=str(group).strip() or None if group else None,
	)


def build_sample_from_record(sample_record: dict) -> Sample:
	"""Create a Sample object from a serialized sample record."""
	if not isinstance(sample_record, dict):
		raise TypeError("sample_record must be a dict")
	return sample_record_to_sample(sample_record)


def samples_from_store(samples_data: Optional[dict]) -> list[Sample]:
	raw = samples_data or {}
	if isinstance(raw, list):
		records = raw
	elif isinstance(raw, dict) and isinstance(raw.get("samples"), list):
		records = raw["samples"]
	elif isinstance(raw, dict):
		records = list(raw.values())
	else:
		records = []
	result: list[Sample] = []
	for record in records:
		if not isinstance(record, dict):
			continue
		try:
			result.append(sample_record_to_sample(record))
		except (KeyError, TypeError):
			continue
	return result


def default_branch(samples: Iterable[Sample], override: Optional[str] = None) -> AnalysisBranch:
	if override:
		return enum_from_value(AnalysisBranch, override, AnalysisBranch.VSEARCH)
	return AnalysisBranch.VSEARCH


def build_run_config(project_data: Optional[dict], samples_data: Optional[dict]) -> RunConfig:
	project_data = resolve_project_data(project_data)
	all_samples = samples_from_store(samples_data)
	samples = [sample for sample in all_samples if getattr(sample, "selected", True)]
	run_id = timestamp_slug()
	branch = default_branch(samples, override=project_data.get("branch"))
	output_dir = project_run_dir(project_data, run_id=run_id)
	analysis_method = str(project_data.get("analysis_method") or "vsearch")
	return RunConfig(
		run_id=str(run_id),
		samples=samples,
		branch=branch,
		analysis_method=analysis_method,
		phylo_mode=PhyloMode.HYBRID,
		output_dir=output_dir,
	)


def run_config_from_store(project_data: Optional[dict], samples_data: Optional[dict], run_config_data: Optional[dict]) -> RunConfig:
	run_config = build_run_config(project_data, samples_data)
	raw = run_config_data or {}
	run_config.enabled_stages = dict(raw.get("enabled_stages", {})) if isinstance(raw.get("enabled_stages"), dict) else {}
	stored_run_id = raw.get("run_id")
	stored_output_dir = raw.get("output_dir")
	if stored_run_id and stored_output_dir:
		run_config.run_id = str(stored_run_id)
		run_config.output_dir = Path(stored_output_dir)
	if raw.get("branch"):
		run_config.branch = enum_from_value(AnalysisBranch, raw.get("branch"), run_config.branch)
	if raw.get("phylo_mode"):
		run_config.phylo_mode = enum_from_value(PhyloMode, raw.get("phylo_mode"), run_config.phylo_mode)
	for key in ("pear_params", "trimmomatic_params", "vsearch_params", "dada2_params", "mafft_params", "phylo_params", "bold_params"):
		value = raw.get(key)
		if isinstance(value, dict):
			current = getattr(run_config, key)
			current.update(value)
	if raw.get("analysis_method"):
		run_config.analysis_method = raw["analysis_method"]
	if raw.get("enabled_stages"):
		run_config.enabled_stages = raw["enabled_stages"]
	method = getattr(run_config, "analysis_method", "vsearch")
	if method == "dada2":
		run_config.branch = AnalysisBranch.DADA2
	elif method == "both":
		run_config.branch = AnalysisBranch.VSEARCH
	else:
		run_config.branch = AnalysisBranch.VSEARCH
	return run_config


def resolve_default_dada2_reference() -> Optional[str]:
	"""Pick a DADA2 assignTaxonomy()-compatible reference fasta.

	PR2's dedicated "dada2" export is the only format actually compatible with
	assignTaxonomy()'s rank-by-';'-split parsing. The raw ARB-SILVA release
	fasta and PR2's mothur fasta use different header conventions and would
	silently misalign taxonomic ranks if passed to assignTaxonomy() directly,
	so they're only used as a last-resort fallback.
	"""
	db_paths = get_db_paths()
	flags = get_feature_flags()
	if flags.enable_pr2 and getattr(db_paths, "pr2_dada2", None) and db_paths.pr2_dada2.exists():
		return str(db_paths.pr2_dada2)
	if flags.enable_pr2 and db_paths.pr2 and db_paths.pr2.exists():
		return str(db_paths.pr2)
	if flags.enable_silva and db_paths.silva and db_paths.silva.exists():
		return str(db_paths.silva)
	return None


def resolve_default_sintax_database(marker: str = "") -> Optional[str]:
	marker_upper = str(marker or "").upper()
	db_paths = get_db_paths()

	# COI: never fall back to SILVA/PR2.
	# Only use a COI-specific SINTAX database if one exists.
	if "COI" in marker_upper or "CO1" in marker_upper:
		for attr in ("coi", "coi_sintax", "bold_sintax"):
			coi_path = getattr(db_paths, attr, None)
			if coi_path and Path(str(coi_path)).exists():
				return str(coi_path)
		return None

	# Non-COI markers can use the usual candidates.
	if "ITS" in marker_upper:
		candidates = [getattr(db_paths, "pr2", None), getattr(db_paths, "silva", None)]
	elif "18S" in marker_upper or "SSU" in marker_upper:
		candidates = [getattr(db_paths, "pr2", None), getattr(db_paths, "silva", None)]
	else:
		candidates = [getattr(db_paths, "pr2", None), getattr(db_paths, "silva", None)]

	for path in candidates:
		if path and Path(str(path)).exists():
			return str(path)

	return None


def serialize_value(value: Any) -> Any:
	if value is None:
		return None
	if isinstance(value, Path):
		return str(value)
	if isinstance(value, dict):
		return {str(key): serialize_value(item) for key, item in value.items()}
	if isinstance(value, (list, tuple, set)):
		return [serialize_value(item) for item in value]
	if hasattr(value, "value") and not is_dataclass(value):
		return value.value
	if is_dataclass(value):
		return serialize_value(asdict(value))
	return value


def serialize_step_output(step_output: StepOutput) -> dict:
	return {
		"step_name": step_output.step_name,
		"success": step_output.success,
		"files": serialize_value(step_output.files),
		"log": serialize_value(step_output.log),
		"metrics": serialize_value(step_output.metrics),
		"error_message": step_output.error_message,
	}


def serialize_run_config(run_config: RunConfig) -> dict:
	return {
		"run_id": run_config.run_id,
		"samples": serialize_value(run_config.samples),
		"branch": run_config.branch.value,
		"analysis_method": getattr(run_config, "analysis_method", "vsearch"),
		"enabled_stages": getattr(run_config, "enabled_stages", {}),
		"phylo_mode": run_config.phylo_mode.value if run_config.phylo_mode else None,
		"output_dir": str(run_config.output_dir),
		"pear_params": serialize_value(run_config.pear_params),
		"trimmomatic_params": serialize_value(run_config.trimmomatic_params),
		"vsearch_params": serialize_value(run_config.vsearch_params),
		"dada2_params": serialize_value(run_config.dada2_params),
		"mafft_params": serialize_value(run_config.mafft_params),
		"phylo_params": serialize_value(run_config.phylo_params),
		"bold_params": serialize_value(run_config.bold_params),
	}


def build_sample_store_entry(sample: Sample, r1_path: Path, r2_path: Optional[Path] = None) -> dict:
	return {
		"internal_id": sample.internal_id,
		"internal_sample_id": sample.internal_id,
		"display_label": sample.display_label,
		"display_name": sample.display_label,
		"mode": sample.mode.value,
		"kind": sample.mode.value,
		"marker": sample.marker.value,
		"r1_path": str(r1_path),
		"r2_path": str(r2_path) if r2_path else None,
		"paths": {
			"r1": str(r1_path),
			"r2": str(r2_path) if r2_path else None,
		},
		"status": "ready",
		"status_color": "green" if sample.mode == AnalysisMode.PE else "teal",
		"selected": sample.selected,
		"is_control": sample.is_control,
		"group": sample.group,
	}
