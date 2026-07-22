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
        "taxonomy_summary_csv": "taxonomy_summary.csv",
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


# Within a phylo/<set>/ directory each artifact lives under its RUNNER's subdir — the
# tree is never at ``phylo/<set>/tree.nwk``. Readers that guessed that flat name found
# nothing after a real run (the pipeline writes ``phylo/<set>/raxml/…``), so the tree,
# the alignment QC and the tip mapping silently blanked. Resolve through the same
# per-set-then-flat ladder the taxonomy files use.
_PHYLO_TREE_CANDIDATES = (
    "raxml/raxml.raxml.support",   # bootstrap support — the viewer colours by this
    "raxml/raxml.raxml.bestTree",  # ML tree, support stripped
    "fasttree/tree.nwk",           # FAST mode, or the HYBRID seed before RAxML
)
_PHYLO_ARTIFACTS = {
    "aligned_fasta": "mafft/aligned.fasta",
    "reversed_asvs_csv": "mafft/reversed_asvs.csv",
}


def phylo_dir(run_output_dir, feature_set: str) -> Optional[Path]:
    """The phylo directory for a feature set: the per-set layout if present, else flat.

    Returns the base ``phylo/`` when the per-set dir is absent so runs that predate
    feature sets (``phylo/raxml/…``) still resolve.
    """
    if not run_output_dir:
        return None
    base = Path(run_output_dir) / "phylo"
    per_set = base / feature_set
    return per_set if per_set.exists() else base


def phylo_tree_path(run_output_dir, feature_set: str) -> Optional[str]:
    """Best available tree for a feature set: bootstrap support, else ML, else FastTree.

    Walks the per-feature-set layout first, then the flat pre-feature-set one, and knows
    the runner subdir each tree sits in — so it returns a path that actually exists.
    """
    if not run_output_dir:
        return None
    base = Path(run_output_dir) / "phylo"
    for root in (base / feature_set, base):
        for rel in _PHYLO_TREE_CANDIDATES:
            candidate = root / rel
            if candidate.exists():
                return str(candidate)
    return None


def phylo_artifact(run_output_dir, feature_set: str, name: str) -> Optional[str]:
    """Path to a named phylo artifact (``aligned_fasta``, ``reversed_asvs_csv``), per-set then flat."""
    rel = _PHYLO_ARTIFACTS.get(name)
    if not rel or not run_output_dir:
        return None
    base = Path(run_output_dir) / "phylo"
    for root in (base / feature_set, base):
        candidate = root / rel
        if candidate.exists():
            return str(candidate)
    return None
