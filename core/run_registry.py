"""The run's own artifacts on disk, as the single source of truth for the UI.

Every selector in Results used to be driven by ``dcc.Store`` data that only the progress
poller ever wrote — and that poller disables itself the moment a run finishes. The stores
are ``storage_type="session"``, so the dropdowns went dead on an app restart, on a new
tab, or simply on arriving at Results without having watched the run go by. The run
directory, meanwhile, has always held the whole truth in ``progress.json``.

Reading from disk makes the UI stateless: it survives F5, navigation, and a restart of
the server, because the facts live where the pipeline put them.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from core.config import get_runs_dir

# Stage key -> the record key the pipeline writes into progress.json.
STAGE_RECORDS = {
    "qc": "qc_record",
    "preprocess": "preprocess_record",
    "analysis": "analysis_record",
    "databases": "db_record",
    "phylo": "phylo_record",
}

STAGE_LABELS = {
    "qc": "QC",
    "preprocess": "Preprocess",
    "analysis": "Analysis",
    "databases": "Databases",
    "phylo": "Phylo",
}

# Appended to a stage's name in a dropdown, so a failed stage is visibly failed rather
# than silently missing from the list.
STATUS_SUFFIX = {
    "success": "",
    "failure": "  ✕ falhou",
    "partial": "  ▲ parcial",
    "skipped": "  — não corrido",
    "pending": "  — pendente",
    "running": "  … a correr",
}

# Stages whose outputs are keyed by sample. Databases and Phylo classify and build trees
# on the POOLED feature set — one OTU table across all samples — so "per sample" is not
# something they have. The selector must say that, not offer an empty list.
PER_SAMPLE_STAGES = frozenset({"qc", "preprocess", "analysis"})


def run_dir_for(run_config_data: Optional[dict] = None, project_data: Optional[dict] = None) -> Optional[Path]:
    """The active run directory: what the UI says, else the newest run on disk.

    The disk fallback is what makes Results work after a restart, when the session
    stores are empty but the run is sitting right there.
    """
    output_dir = (run_config_data or {}).get("output_dir")
    if output_dir and (Path(output_dir) / "progress.json").exists():
        return Path(output_dir)
    return latest_run_dir(project_data)


def latest_run_dir(project_data: Optional[dict] = None) -> Optional[Path]:
    """Newest run directory that actually has a progress.json, or ``None``."""
    runs_root = Path(get_runs_dir())
    if not runs_root.exists():
        return None

    candidates: list[Path] = []
    # Runs live at <runs>/<project>/<timestamp>/. Prefer the active project's, but fall
    # back to every project so a user who never opened Setup still sees their last run.
    search_roots = []
    active = (project_data or {}).get("active_id")
    projects = (project_data or {}).get("projects") or {}
    name = (projects.get(active) or {}).get("project_name") if active else None
    if name:
        search_roots.append(runs_root / str(name))
    search_roots.append(runs_root)

    for root in search_roots:
        if not root.exists():
            continue
        candidates = [path.parent for path in root.glob("*/progress.json")]
        if not candidates:
            candidates = [path.parent for path in root.glob("*/*/progress.json")]
        if candidates:
            break

    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def load_run_state(run_dir) -> dict:
    """The run's progress.json, or ``{}``."""
    if not run_dir:
        return {}
    path = Path(run_dir) / "progress.json"
    if not path.exists():
        return {}
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return state if isinstance(state, dict) else {}


def stage_record(state: dict, stage: str) -> dict:
    record = (state or {}).get(STAGE_RECORDS.get(stage, ""), None)
    return record if isinstance(record, dict) else {}


def run_stages(state: dict) -> list[dict]:
    """Every stage with its real status — including the ones that failed or were skipped.

    Listing only successful stages (which is what Results did) is why the dropdown was
    empty on the run the user actually cared about: analysis and databases had failed,
    so there was nothing to select and no explanation of why.
    """
    stages = []
    for key, label in STAGE_LABELS.items():
        record = stage_record(state, key)
        status = str(record.get("status") or "pending")
        stages.append(
            {
                "key": key,
                "label": label,
                "status": status,
                "has_outputs": bool(record.get("outputs")),
                "per_sample": key in PER_SAMPLE_STAGES,
            }
        )
    return stages


def list_run_samples(state: dict, stage: Optional[str] = None) -> list[dict]:
    """``[{"label": ..., "value": ...}, ...]`` for the samples this run produced.

    Read from the run's own outputs, not from the samples store: the store reflects what
    the user last configured, which is not necessarily what the run actually processed.
    """
    seen: dict[str, str] = {}

    stages = [stage] if stage else list(STAGE_LABELS)
    for key in stages:
        outputs = stage_record(state, key).get("outputs") or {}
        for sample_id, value in outputs.items():
            # Stage-level artifacts (combined_fasta, taxonomy, feature_sets, ...) live in
            # the same dict as the per-sample records; only the latter are dicts of tools
            # or StepOutputs keyed by a sample id.
            if not isinstance(value, dict) or sample_id in _NON_SAMPLE_KEYS:
                continue
            seen.setdefault(sample_id, sample_id)

    # The preprocess record carries the display labels the user recognises.
    for entry in stage_record(state, "preprocess").get("next_samples") or []:
        if isinstance(entry, dict) and entry.get("internal_id") in seen:
            seen[entry["internal_id"]] = entry.get("display_label") or entry["internal_id"]

    return [{"label": label, "value": sample_id} for sample_id, label in sorted(seen.items())]


_NON_SAMPLE_KEYS = frozenset(
    {
        "combined_fasta",
        "combined_fasta_vsearch",
        "combined_fasta_dada2",
        "dada2",
        "taxonomy",
        "taxonomy_by_set",
        "phylo",
        "phylo_by_set",
        "vsearch_pool",
        "feature_sets",
        "feature_set_comparison",
        "branch_errors",
        "dropped_samples",
        "otu_table_tsv",
        "otu_matrix",
        "combined_otu_uc",
    }
)


def overall_status(state: dict) -> str:
    """One badge for the whole run: success / partial / failure / running."""
    if not state:
        return "pending"
    statuses = [str(stage_record(state, key).get("status") or "pending") for key in STAGE_LABELS]
    ran = [status for status in statuses if status not in ("pending", "skipped")]
    if not ran:
        return "pending"
    if not state.get("done"):
        return "running"
    if all(status == "success" for status in ran):
        return "success"
    if any(status == "success" for status in ran):
        # Half the pipeline delivering is a materially different outcome from none of it,
        # and the user needs to be told which one they have.
        return "partial"
    return "failure"
