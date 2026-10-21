"""Persisted phylogenetics metrics: alignment QC, branch support, UniFrac.

Computed once in the Phylo stage. UniFrac in particular is O(tips x samples) per
pair and must not be redone on every page render.

Branch support means different things depending on the tree builder: FastTree
writes local support on 0-1, RAxML-NG writes bootstrap on 0-100. Both are
normalised to a percentage here so a single 70% threshold reads the same either
way.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Optional

import numpy as np

from analysis.analysis_metrics import load_matrix
from phylo.phylo_diversity import load_asv_mapping

# The conventional cut-off below which a clade is not considered supported.
SUPPORT_THRESHOLD = 70.0

# MAFFT --adjustdirection prefixes reverse-complemented records with this.
_REVERSED_PREFIX = "_R_"

_BASES = ("A", "C", "G", "T")


def read_alignment(aligned_fasta) -> tuple[list[str], list[str]]:
    """Return ``(ids, sequences)`` from an aligned FASTA, sequences upper-cased."""
    path = Path(aligned_fasta)
    if not path.exists() or path.stat().st_size == 0:
        return [], []
    ids: list[str] = []
    sequences: list[str] = []
    current: list[str] = []
    with path.open("r", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
                ids.append(line[1:].strip())
            else:
                current.append(line)
    if current:
        sequences.append("".join(current).upper())
    return ids, sequences


def reversed_sequences(ids: list[str]) -> list[str]:
    """Ids MAFFT reverse-complemented. Many of them means mixed-orientation input."""
    return [name for name in ids if name.startswith(_REVERSED_PREFIX)]


def alignment_summary(aligned_fasta) -> dict:
    """Sequence count, alignment length, overall gap share and reversed-sequence count."""
    ids, sequences = read_alignment(aligned_fasta)
    if not sequences:
        return {}
    length = max(len(sequence) for sequence in sequences)
    total = length * len(sequences)
    gaps = sum(sequence.count("-") for sequence in sequences)
    reversed_ids = reversed_sequences(ids)
    return {
        "n_sequences": len(sequences),
        "alignment_length": length,
        "pct_gaps": round(100.0 * gaps / total, 3) if total else 0.0,
        "n_reversed": len(reversed_ids),
        "reversed_ids": ";".join(reversed_ids),
    }


def alignment_columns(aligned_fasta) -> list[dict]:
    """Per-column gap fraction and conservation.

    Conservation is the frequency of the most common *non-gap* base, so a column
    that is 90% gaps and 10% a single base scores 1.0 on conservation and 0.9 on
    gaps — the two numbers have to be read together.
    """
    _, sequences = read_alignment(aligned_fasta)
    if not sequences:
        return []
    length = max(len(sequence) for sequence in sequences)
    depth = len(sequences)
    rows: list[dict] = []
    for position in range(length):
        column = [sequence[position] if position < len(sequence) else "-" for sequence in sequences]
        gaps = column.count("-")
        residues = Counter(base for base in column if base in _BASES)
        occupancy = sum(residues.values())
        conservation = (residues.most_common(1)[0][1] / occupancy) if occupancy else 0.0
        rows.append(
            {
                "position": position + 1,
                "gap_fraction": round(gaps / depth, 4),
                "conservation": round(conservation, 4),
                "occupancy": occupancy,
            }
        )
    return rows


# ── Branch support ───────────────────────────────────────────────────────


def branch_support(tree_path) -> list[dict]:
    """Support value of every internal node, normalised to a percentage.

    FastTree writes local support on 0-1 and RAxML-NG bootstrap on 0-100. Scaling
    is decided from the whole tree, not per node: a single 0.86 among values that
    otherwise run to 100 would otherwise be rescaled to 86%.
    """
    from Bio import Phylo

    path = Path(tree_path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    try:
        tree = Phylo.read(str(path), "newick")
    except Exception:
        return []

    raw: list[float] = []
    for clade in tree.get_nonterminals():
        value = clade.confidence
        if value is None and clade.name:
            try:
                value = float(clade.name)
            except ValueError:
                value = None
        if value is not None:
            raw.append(float(value))
    if not raw:
        return []
    scale = 100.0 if max(raw) <= 1.0 else 1.0
    return [
        {"node": index, "support": round(value * scale, 2), "is_supported": int(value * scale >= SUPPORT_THRESHOLD)}
        for index, value in enumerate(raw, start=1)
    ]


# ── UniFrac ──────────────────────────────────────────────────────────────


def unifrac_distances(tree_path, matrix_path, id_mapping_csv=None, weighted: bool = True) -> dict:
    """UniFrac distances between samples, keyed by the abundance table's samples.

    Feature ids in the matrix are the original OTU ids; the tree's tips were
    renamed to ``ASV1..N`` before alignment, so the mapping is mandatory —
    without it no tip matches any feature and every distance is zero.
    """
    from skbio import TreeNode
    from skbio.diversity import beta_diversity

    tree_path = Path(tree_path)
    if not tree_path.exists() or not (matrix_path and Path(matrix_path).exists()):
        return {}
    features, samples, counts = load_matrix(matrix_path)
    if not features or len(samples) < 2:
        return {}

    mapping = load_asv_mapping(id_mapping_csv) if id_mapping_csv else {}
    if not mapping:
        return {}

    try:
        tree = TreeNode.read(str(tree_path))
        # UniFrac needs a rooted tree; FastTree emits an unrooted trifurcation.
        tree = tree.root_at_midpoint()
    except Exception:
        return {}

    # Keep only features whose tip actually exists in this tree. scikit-bio raises
    # on a taxon it cannot find, and skipping the check would let a stale mapping
    # produce a matrix of silent zeros.
    tips = {node.name for node in tree.tips() if node.name}
    taxa = [mapping.get(feature.split(";")[0].strip()) for feature in features]
    # Reverse-complemented records keep MAFFT's "_R_" prefix on the tree's tips.
    taxa = [name if (not name or name in tips) else (_REVERSED_PREFIX + name) for name in taxa]
    keep = [index for index, name in enumerate(taxa) if name and name in tips]
    if not keep:
        return {}

    try:
        metric = "weighted_unifrac" if weighted else "unweighted_unifrac"
        distances = beta_diversity(
            metric,
            counts[:, keep].astype(int),
            ids=samples,
            taxa=[taxa[index] for index in keep],
            tree=tree,
        )
    except Exception:
        return {}
    return {"samples": list(samples), "matrix": np.asarray(distances.data).round(6).tolist()}


# ── Persistence ──────────────────────────────────────────────────────────


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_distance_csv(path: Path, payload: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", *payload["samples"]])
        for name, row in zip(payload["samples"], payload["matrix"]):
            writer.writerow([name, *row])
    return path


def build_phylo_metrics(
    phylo_dir,
    tree_path=None,
    aligned_fasta=None,
    matrix_path=None,
    id_mapping_csv=None,
) -> dict[str, Path]:
    """Compute and persist the phylo-stage tables. Returns name -> path."""
    phylo_dir = Path(phylo_dir)
    metrics_dir = phylo_dir / "metrics"
    written: dict[str, Path] = {}

    if aligned_fasta and Path(aligned_fasta).exists():
        summary = alignment_summary(aligned_fasta)
        if summary:
            written["alignment_summary_csv"] = _write_csv(metrics_dir / "alignment_summary.csv", list(summary), [summary])
        columns = alignment_columns(aligned_fasta)
        if columns:
            written["alignment_columns_csv"] = _write_csv(
                metrics_dir / "alignment_columns.csv", ["position", "gap_fraction", "conservation", "occupancy"], columns
            )

    if tree_path and Path(tree_path).exists():
        support = branch_support(tree_path)
        if support:
            written["branch_support_csv"] = _write_csv(
                metrics_dir / "branch_support.csv", ["node", "support", "is_supported"], support
            )

        for weighted, key, name in ((True, "weighted_unifrac_csv", "weighted_unifrac.csv"), (False, "unweighted_unifrac_csv", "unweighted_unifrac.csv")):
            payload = unifrac_distances(tree_path, matrix_path, id_mapping_csv, weighted=weighted)
            if payload:
                written[key] = _write_distance_csv(metrics_dir / name, payload)

    return written


def support_tree_path(phylo_files: dict, phylo_dir=None) -> Optional[Path]:
    """Pick the tree that actually carries branch support.

    RAxML-NG writes support values to ``*.support`` and leaves ``bestTree`` bare,
    but ``tree_nwk`` is aliased to ``bestTree`` for the download link. Branch
    support and the support-coloured tree must therefore read ``support_tree``
    first; FastTree runs have no separate support file, so ``tree_nwk`` is right
    there.
    """
    # The support tree wins first: bestTree/tree_nwk are bare, so falling back to a
    # real ``*.support`` file on disk must take priority over them. The record's
    # ``support_tree`` can be stale (e.g. pointing at ``*.supportFBP`` when RAxML
    # wrote ``*.support``), which is why the on-disk probe is not merely a last
    # resort.
    recorded_support = phylo_files.get("support_tree")
    if recorded_support and Path(str(recorded_support)).exists():
        return Path(str(recorded_support))
    if phylo_dir is not None:
        for name in ("raxml.raxml.support", "raxml.raxml.supportFBP", "raxml.raxml.supportTBE"):
            candidate = Path(phylo_dir) / "raxml" / name
            if candidate.exists():
                return candidate
    for value in (phylo_files.get("tree_nwk"), phylo_files.get("best_tree")):
        if value and Path(str(value)).exists():
            return Path(str(value))
    if phylo_dir is not None:
        candidate = Path(phylo_dir) / "fasttree" / "tree.nwk"
        if candidate.exists():
            return candidate
    return None


def build_metrics_from_files(phylo_dir, phylo_files: dict, matrix_path=None) -> dict[str, str]:
    """Persist the phylo tables from a stage's file map.

    Shared by the pipeline thread and the page's own "Run Phylo" button so both
    write the same tables from the same artifacts.
    """
    phylo_dir = Path(phylo_dir)

    def _path(value, fallback: str):
        if value and Path(str(value)).exists():
            return Path(str(value))
        candidate = phylo_dir / fallback
        return candidate if candidate.exists() else None

    written = build_phylo_metrics(
        phylo_dir,
        tree_path=support_tree_path(phylo_files, phylo_dir),
        aligned_fasta=_path(phylo_files.get("aligned_fasta"), "mafft/aligned.fasta"),
        matrix_path=matrix_path,
        id_mapping_csv=_path(phylo_files.get("reversed_asvs_csv"), "mafft/reversed_asvs.csv"),
    )
    return {key: str(path) for key, path in written.items()}


def phylo_metrics_paths(phylo_dir) -> dict[str, Path]:
    """Return the persisted phylo tables that exist under ``phylo_dir``."""
    metrics_dir = Path(phylo_dir) / "metrics"
    names = {
        "alignment_summary_csv": "alignment_summary.csv",
        "alignment_columns_csv": "alignment_columns.csv",
        "branch_support_csv": "branch_support.csv",
        "weighted_unifrac_csv": "weighted_unifrac.csv",
        "unweighted_unifrac_csv": "unweighted_unifrac.csv",
    }
    return {key: metrics_dir / name for key, name in names.items() if (metrics_dir / name).exists()}
