"""Storage, retention, cleanup and session metadata helpers."""

from __future__ import annotations

import gzip
import os
import shutil
import time
import zipfile
from datetime import datetime

from .config import DELETE_AFTER_STEPS, KEEP_INTERMEDIATES
from .paths import PROJECTS_ROOT, _project_exports_dir, _project_session_path, _run_session_path
from .utils import _path_is_valid, _read_json, _write_json


class StorageConfig:
    """Configuration for storage management policies."""

    MAX_RUNS_TO_KEEP = int(os.getenv("MAX_RUNS_TO_KEEP", "1"))
    MAX_RUN_AGE_DAYS = int(os.getenv("MAX_RUN_AGE_DAYS", "30"))
    AUTO_DELETE_TEMP = str(os.getenv("AUTO_DELETE_TEMP", "1")).strip().lower() in {"1", "true", "yes", "on"}
    TEMP_FILE_MAX_AGE_DAYS = int(os.getenv("TEMP_FILE_MAX_AGE_DAYS", "3"))
    DELETE_INTERMEDIATE_STEPS = [
        step.strip()
        for step in os.getenv("DELETE_INTERMEDIATE_STEPS", "qc,preprocess").split(",")
        if step.strip()
    ]
    MAX_EXPORTS_TO_KEEP = int(os.getenv("MAX_EXPORTS_TO_KEEP", "10"))


def _step_manifest_path(project_state, step_name):
    from .paths import _resolve_output_path

    return _resolve_output_path(project_state, step_name, "final", f"{step_name}_manifest.json")


def _final_manifest_path(project_state):
    from .paths import _resolve_output_path

    return _resolve_output_path(project_state, "analysis", "final", "final_manifest.json")


def _step_output_entry(path, label=None):
    return {
        "label": label or os.path.basename(path) if path else "",
        "path": path or "",
        "exists": bool(path) and os.path.exists(path),
        "size": int(os.path.getsize(path)) if path and os.path.exists(path) else 0,
        "valid": _path_is_valid(path),
    }


def _build_step_manifest(step_name, status, outputs=None, inputs=None, error=None, notes=None):
    return {
        "step": step_name,
        "status": status,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "outputs": outputs or {},
        "inputs": inputs or {},
        "error": error or "",
        "notes": notes or "",
    }


def _write_step_manifest(project_state, step_name, manifest):
    if not isinstance(manifest, dict):
        return False
    return _write_json(_step_manifest_path(project_state, step_name), manifest)


def _read_step_manifest(project_state, step_name, default=None):
    return _read_json(_step_manifest_path(project_state, step_name), default=default)


def _manifest_path_value(payload, *keys):
    current = payload or {}
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current if _path_is_valid(current) else None


def _get_all_project_dirs():
    projects = []
    if not os.path.exists(PROJECTS_ROOT):
        return projects
    try:
        for entry in os.listdir(PROJECTS_ROOT):
            pdir = os.path.join(PROJECTS_ROOT, entry)
            if os.path.isdir(pdir):
                projects.append(pdir)
    except Exception:
        return projects
    return projects


def _latest_run_dir(project_dir):
    runs_dir = os.path.join(project_dir, "runs")
    if not os.path.exists(runs_dir):
        return None
    candidates = []
    for entry in os.listdir(runs_dir):
        rdir = os.path.join(runs_dir, entry)
        if os.path.isdir(rdir):
            candidates.append((os.path.getmtime(rdir), rdir))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates[0][1]


def _build_run_session_metadata(project_state, run_config, workflow_scope, samples_store, completed_at=None, outputs=None):
    state = project_state or {}
    workflow_samples = _get_workflow_scope_samples(workflow_scope)
    if not workflow_samples:
        workflow_samples = list((samples_store or {}).keys())

    return {
        "run_id": os.path.basename(state.get("run_dir") or ""),
        "project_name": state.get("project_name", ""),
        "project_dir": state.get("project_dir", ""),
        "run_dir": state.get("run_dir", ""),
        "created_at": (run_config or {}).get("created_at") or datetime.now().isoformat(timespec="seconds"),
        "completed_at": completed_at or "",
        "samples_selected": workflow_samples,
        "steps_executed": (run_config or {}).get("steps", []),
        "parameters": run_config or {},
        "tool_versions": _collect_tool_versions(),
        "outputs": outputs or {},
        "execution_log": [],
    }


def _write_run_session_metadata(project_state, run_config, workflow_scope, samples_store, completed_at=None, outputs=None):
    state = project_state or {}
    run_dir = state.get("run_dir")
    if not run_dir:
        return False
    metadata = _build_run_session_metadata(state, run_config, workflow_scope, samples_store, completed_at=completed_at, outputs=outputs)
    try:
        _write_json(_run_session_path(run_dir), metadata)
        return True
    except Exception:
        return False


def _update_run_session_metadata(project_state, updates):
    state = project_state or {}
    run_dir = state.get("run_dir")
    if not run_dir:
        return False
    session_json_path = _run_session_path(run_dir)
    session_payload = _read_json(session_json_path, default={}) or {}
    if not isinstance(session_payload, dict):
        session_payload = {}
    if not isinstance(updates, dict) or not updates:
        return False
    session_payload.update(updates)
    try:
        _write_json(session_json_path, session_payload)
        return True
    except Exception:
        return False


def _list_existing_projects():
    projects = []
    if not os.path.exists(PROJECTS_ROOT):
        return projects
    for entry in sorted(os.listdir(PROJECTS_ROOT)):
        pdir = os.path.join(PROJECTS_ROOT, entry)
        if os.path.isdir(pdir) and os.path.exists(os.path.join(pdir, "project.json")):
            projects.append({"label": entry, "value": entry})
    return projects


def _apply_cleanup_policy(run_dir, delete_after_steps=None):
    if KEEP_INTERMEDIATES or not run_dir:
        return []
    steps_to_delete = [step for step in (delete_after_steps if delete_after_steps is not None else DELETE_AFTER_STEPS) if step]
    deleted = []
    for step_name in steps_to_delete:
        step_dir = os.path.join(run_dir, "steps", step_name)
        if os.path.exists(step_dir):
            shutil.rmtree(step_dir, ignore_errors=True)
            deleted.append(step_name)
    return deleted


def _apply_cleanup_post_run(project_dir, run_dir):
    """Apply cleanup after a run completes."""
    if not project_dir or not isinstance(project_dir, str):
        return False
    if not run_dir or not isinstance(run_dir, str):
        return False
    try:
        StorageManager.cleanup_old_runs(project_dir)
        if StorageConfig.AUTO_DELETE_TEMP:
            StorageManager.cleanup_temp_files(project_dir)
        StorageManager.compact_run_outputs(run_dir)
        StorageManager.cleanup_export_archives(project_dir, run_dir)
        return True
    except Exception:
        return False


class StorageManager:
    """Manages storage cleanup and retention policies."""

    @staticmethod
    def cleanup_old_runs(project_dir):
        runs_dir = os.path.join(project_dir, "runs")
        if not os.path.exists(runs_dir):
            return []
        deleted = []
        try:
            entries = []
            for entry in os.listdir(runs_dir):
                entry_path = os.path.join(runs_dir, entry)
                if os.path.isdir(entry_path):
                    entries.append((os.path.getmtime(entry_path), entry_path, entry))
            entries.sort(reverse=True)
            for _, entry_path, entry in entries[StorageConfig.MAX_RUNS_TO_KEEP:]:
                try:
                    shutil.rmtree(entry_path)
                    deleted.append(entry)
                except Exception:
                    pass
        except Exception:
            pass
        return deleted

    @staticmethod
    def cleanup_temp_files(project_dir):
        deleted = []
        temp_patterns = [
            ("steps/*/temp_*", "step temp files"),
            ("runs/*/steps/*/tmp_*", "run temp files"),
            ("uploads/*_temp_*", "upload temp files"),
        ]
        now = time.time()
        max_age_seconds = StorageConfig.TEMP_FILE_MAX_AGE_DAYS * 86400
        try:
            import glob

            for pattern, _ in temp_patterns:
                full_pattern = os.path.join(project_dir, pattern)
                for filepath in glob.glob(full_pattern):
                    try:
                        if now - os.path.getmtime(filepath) > max_age_seconds:
                            if os.path.isfile(filepath):
                                os.remove(filepath)
                            elif os.path.isdir(filepath):
                                shutil.rmtree(filepath)
                            deleted.append(os.path.basename(filepath))
                    except Exception:
                        pass
        except Exception:
            pass
        return deleted

    @staticmethod
    def compact_run_outputs(run_dir):
        compressed_files = []
        try:
            import glob

            for pattern in ("steps/*/intermediate/*.fastq", "steps/*/intermediate/*.fq"):
                full_pattern = os.path.join(run_dir, pattern)
                for filepath in glob.glob(full_pattern):
                    try:
                        if os.path.getsize(filepath) > 10 * 1024 * 1024:
                            gz_path = filepath + ".gz"
                            if not os.path.exists(gz_path):
                                with open(filepath, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                                os.remove(filepath)
                                compressed_files.append(os.path.basename(filepath))
                    except Exception:
                        pass
        except Exception:
            pass
        return compressed_files

    @staticmethod
    def cleanup_export_archives(project_dir, run_dir=None):
        deleted = []
        try:
            if run_dir and os.path.isdir(run_dir):
                export_staging_dir = os.path.join(run_dir, "export_staging")
                if os.path.isdir(export_staging_dir):
                    shutil.rmtree(export_staging_dir, ignore_errors=True)
                for name in os.listdir(run_dir):
                    file_path = os.path.join(run_dir, name)
                    if os.path.isfile(file_path) and name.lower().endswith(".zip"):
                        try:
                            os.remove(file_path)
                            deleted.append(file_path)
                        except OSError:
                            pass

            exports_dir = _project_exports_dir(project_dir) if project_dir else None
            if exports_dir and os.path.isdir(exports_dir):
                zip_entries = []
                for name in os.listdir(exports_dir):
                    file_path = os.path.join(exports_dir, name)
                    if os.path.isfile(file_path) and name.lower().endswith(".zip"):
                        try:
                            zip_entries.append((os.path.getmtime(file_path), file_path))
                        except OSError:
                            continue

                zip_entries.sort(key=lambda item: item[0], reverse=True)
                for _, file_path in zip_entries[StorageConfig.MAX_EXPORTS_TO_KEEP :]:
                    try:
                        os.remove(file_path)
                        deleted.append(file_path)
                    except OSError:
                        pass
        except Exception:
            pass
        return deleted


def _load_project_session(project_dir):
    payload = _read_json(_project_session_path(project_dir), default={}) or {}
    state = payload.get("project_state", {}) if isinstance(payload, dict) else {}
    state = state if isinstance(state, dict) else {}
    run_config = payload.get("run_config", {}) if isinstance(payload, dict) else {}
    run_config = run_config if isinstance(run_config, dict) else {}
    uploaded_files = payload.get("uploaded_files", {"preprocess_mode": "se"}) if isinstance(payload, dict) else {"preprocess_mode": "se"}
    uploaded_files = uploaded_files if isinstance(uploaded_files, dict) else {"preprocess_mode": "se"}
    active_sample = payload.get("active_sample", "__ALL__") if isinstance(payload, dict) else "__ALL__"
    samples_store = payload.get("samples_store", {}) if isinstance(payload, dict) else {}
    samples_store = samples_store if isinstance(samples_store, dict) else {}
    workflow_scope = payload.get("workflow_scope", {}) if isinstance(payload, dict) else {}
    workflow_scope = workflow_scope if isinstance(workflow_scope, dict) else {}
    return state, run_config, uploaded_files, active_sample, samples_store, workflow_scope


def _save_project_session(project_state, run_config, uploaded_files, active_sample, samples_store, workflow_scope):
    state = project_state or {}
    project_dir = state.get("project_dir")
    if not project_dir:
        return False
    payload = {
        "saved_at": datetime.now().isoformat(timespec="seconds"),
        "project_state": state,
        "run_config": run_config or {},
        "uploaded_files": uploaded_files or {"preprocess_mode": "se"},
        "active_sample": active_sample if active_sample is not None else "__ALL__",
        "samples_store": samples_store or {},
        "workflow_scope": workflow_scope or {},
    }
    _write_json(_project_session_path(project_dir), payload)
    return True


def _get_workflow_scope_samples(workflow_scope):
    scope = workflow_scope or {}
    samples = scope.get("samples_selected", [])
    if not samples or "__ALL__" in samples:
        return []
    return samples


def _collect_tool_versions():
    # Technical debt: this currently returns an empty payload.
    # Older workflow implementations queried tool binaries via subprocess
    # for reproducibility metadata; restore that behavior when wiring this
    # module fully into the runtime execution path.
    return {}
