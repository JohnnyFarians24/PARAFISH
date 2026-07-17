"""Server-side persistence for projects and their samples.

Projects live in the browser store (localStorage) for responsiveness, but that
is per-browser and lost on a different machine or after clearing site data. This
module mirrors each project (metadata + its sample store) to
``projects/<id>/project.json`` so previously created projects can be listed and
reopened from any session, with their uploaded samples restored.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


def _projects_root() -> Path:
    root = Path("projects")
    root.mkdir(parents=True, exist_ok=True)
    return root


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def save_project_bundle(project_id: str, project: dict, samples: Optional[dict] = None) -> bool:
    """Persist a project's metadata and sample store to disk."""
    if not project_id:
        return False
    try:
        project_dir = _projects_root() / str(project_id)
        project_dir.mkdir(parents=True, exist_ok=True)
        payload = {"project": project or {}, "samples": samples or {}, "saved_at": _now()}
        (project_dir / "project.json").write_text(json.dumps(payload, default=str), encoding="utf-8")
        return True
    except Exception:
        return False


def load_all_projects() -> dict:
    """Return ``{project_id: project_metadata}`` for every persisted project."""
    projects: dict = {}
    root = _projects_root()
    if not root.exists():
        return projects
    for entry in sorted(root.iterdir()):
        manifest = entry / "project.json"
        if not (entry.is_dir() and manifest.exists()):
            continue
        try:
            payload = json.loads(manifest.read_text(encoding="utf-8"))
        except Exception:
            continue
        project = payload.get("project") or {}
        project_id = project.get("id") or entry.name
        project.setdefault("id", project_id)
        projects[project_id] = project
    return projects


def load_project_bundle(project_id: str) -> tuple[Optional[dict], Optional[dict]]:
    """Return ``(project_metadata, sample_store)`` for a persisted project."""
    manifest = _projects_root() / str(project_id) / "project.json"
    if not manifest.exists():
        return None, None
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except Exception:
        return None, None
    return payload.get("project") or {}, payload.get("samples") or {}


def delete_project(project_id: str) -> bool:
    project_dir = _projects_root() / str(project_id)
    if not project_dir.exists():
        return False
    try:
        shutil.rmtree(project_dir)
        return True
    except Exception:
        return False
