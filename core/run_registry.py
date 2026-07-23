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


def _run_dir_is_complete(run_dir: Path) -> bool:
    state = load_run_state(run_dir)
    return bool(state.get("done")) and not state.get("interrupted")


def newest_complete_run(run_dirs) -> Optional[Path]:
    """The newest FINISHED run among ``run_dirs`` — the SINGLE source of truth for
    "which run does the UI restore", shared by Results (``run_dir_for``) and the Eixo B
    pages (``ui.runtime_state.latest_run_dir``) so the two can never pick different runs
    in one session.

    Ordered by run_id — the immutable UTC directory name, whose lexical order IS
    chronological — and NEVER by mtime, which a re-freeze, a copy, or an editor touching
    a file inside an old run bumps above the real newest. Skips OOM orphans (no ``done``,
    or ``interrupted``): their partial records render the pages blank, which reads as the
    restore having broken.
    """
    complete = [Path(run_dir) for run_dir in run_dirs if _run_dir_is_complete(Path(run_dir))]
    return max(complete, key=lambda path: path.name) if complete else None


def latest_run_dir(project_data: Optional[dict] = None) -> Optional[Path]:
    """Newest COMPLETE run directory for the UI, or ``None``.

    Runs live at <runs>/<project>/<timestamp>/. Prefer the active project's, but fall
    back to every project so a user who never opened Setup still sees their last run —
    the one axis on which this legitimately differs from the project-scoped Eixo B
    fallback. The run SELECTION (newest complete by run_id) is the shared primitive, so
    Results and Eixo B agree on which run within a scope.
    """
    runs_root = Path(get_runs_dir())
    if not runs_root.exists():
        return None

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
        run_dirs = [path.parent for path in root.glob("*/progress.json")]
        if not run_dirs:
            run_dirs = [path.parent for path in root.glob("*/*/progress.json")]
        if run_dirs:
            # This root owns the answer (the first with any runs), so its runs decide it;
            # if they are all orphans the result is None rather than leaking to another
            # project's runs.
            return newest_complete_run(run_dirs)
    return None


def run_claims_running(state: dict) -> bool:
    """A run that has not finished but still shows a stage as ``running``.

    That is the exact fingerprint of an OOM-killed / abandoned run: the process that
    would have written ``done`` and flipped the stage to a terminal status is gone, so
    the badge lies "running" forever.
    """
    if not state or state.get("done") or state.get("interrupted"):
        return False
    stages = state.get("stages") or {}
    return any(status == "running" for status in stages.values())


def reconcile_interrupted_runs(runs_root=None, is_live=None) -> list[Path]:
    """Mark every run that claims to be running but has no live owning thread as
    interrupted, and return the run dirs touched.

    A pipeline runs in a thread inside the server process, so a run left "running" on
    disk is one whose process died (OOM kill or restart) — there is no way for a run to
    outlive its process. This is called at server startup, where no pipeline thread can
    have survived, so ``is_live`` is universally False; the predicate exists so the
    invariant ("only reconcile what has no live thread") is explicit and testable, and
    so a caller that runs it while the server is up cannot clobber a genuinely live run.

    Reconciliation writes ``interrupted: true`` and flips any ``running`` stage to
    ``interrupted``. It deliberately does NOT set ``done: true``: the run did not finish,
    so it must stay out of the completed-run set that drives restore.
    """
    root = Path(runs_root) if runs_root else Path(get_runs_dir())
    if not root.exists():
        return []
    is_live = is_live or (lambda run_id: False)
    reconciled: list[Path] = []
    seen: set[Path] = set()
    # Runs live at <runs>/<project>/<run_id>/; some legacy layouts are one level up.
    for pattern in ("*/progress.json", "*/*/progress.json"):
        for progress_path in root.glob(pattern):
            if progress_path in seen:
                continue
            seen.add(progress_path)
            state = load_run_state(progress_path.parent)
            if not run_claims_running(state):
                continue
            if is_live(state.get("run_id")):
                continue
            state["interrupted"] = True
            state["overall"] = "interrupted"
            state["stages"] = {
                name: ("interrupted" if status == "running" else status)
                for name, status in (state.get("stages") or {}).items()
            }
            progress_path.write_text(json.dumps(state), encoding="utf-8")
            reconciled.append(progress_path.parent)
    return reconciled


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
