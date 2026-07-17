"""Resolving which feature set (OTU or ASV) a page is showing.

Every stage now writes one directory per feature set — ``databases/otu/``,
``databases/asv/``, ``phylo/otu/`` … — so a page must always know WHICH set it is
rendering. Guessing, or falling back to "whatever file exists", is how an ASV count
ends up under an OTU heading.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

FEATURE_SET_LABELS = {"otu": "OTUs (VSEARCH)", "asv": "ASVs (DADA2)"}
FEATURE_SET_NOUN = {"otu": "OTU", "asv": "ASV"}

DEFAULT_FEATURE_SET = "otu"


def available_feature_sets(analysis_outputs: Optional[dict], run_output_dir=None) -> list[dict]:
    """Feature sets this run actually produced, as ``[{"key", "label"}, ...]``."""
    recorded = (analysis_outputs or {}).get("feature_sets") or []
    sets = [
        {"key": entry["key"], "label": entry.get("label") or FEATURE_SET_LABELS.get(entry["key"], entry["key"])}
        for entry in recorded
        if isinstance(entry, dict) and entry.get("key")
    ]
    if sets:
        return sets
    # Runs that predate feature sets, or a page loaded straight from disk.
    if run_output_dir:
        found = []
        for key in ("otu", "asv"):
            if (Path(run_output_dir) / "databases" / key).exists() or (Path(run_output_dir) / "analysis" / "metrics" / key).exists():
                found.append({"key": key, "label": FEATURE_SET_LABELS[key]})
        if found:
            return found
    return [{"key": DEFAULT_FEATURE_SET, "label": FEATURE_SET_LABELS[DEFAULT_FEATURE_SET]}]


def resolve_feature_set(selected: Optional[str], analysis_outputs: Optional[dict], run_output_dir=None) -> str:
    """The feature set to render: the user's choice when it is available, else the first."""
    sets = available_feature_sets(analysis_outputs, run_output_dir)
    keys = [entry["key"] for entry in sets]
    if selected and selected in keys:
        return selected
    return keys[0] if keys else DEFAULT_FEATURE_SET


def taxonomy_file(db_outputs: Optional[dict], run_output_dir, feature_set: str, name: str) -> Optional[str]:
    """Path to one of the Databases stage's per-feature-set outputs, or ``None``.

    ``name`` is a key such as ``taxonomy_consensus_csv`` or ``taxonomy_hits_csv``.
    """
    by_set = ((db_outputs or {}).get("taxonomy_by_set") or {}).get(feature_set) or {}
    value = (by_set.get("files") or {}).get(name)
    if isinstance(value, str) and value and Path(value).exists():
        return value

    if not run_output_dir:
        return None
    filename = {
        "taxonomy_consensus_csv": "taxonomy_consensus.csv",
        "taxonomy_hits_csv": "taxonomy_hits.csv",
        "source_status_csv": "source_status.csv",
        "conflict_report_csv": "conflict_report.csv",
        "source_assignments_csv": "source_assignments.csv",
        "annotated_fasta": "annotated.fasta",
    }.get(name)
    if not filename:
        return None
    # Per-feature-set layout first; then the flat pre-feature-set layout.
    for candidate in (
        Path(run_output_dir) / "databases" / feature_set / filename,
        Path(run_output_dir) / "databases" / filename,
    ):
        if candidate.exists():
            return str(candidate)
    return None
