"""Phylogenetic diversity and tree-layout helpers.

- ``faith_pd``: rooted Faith's Phylogenetic Diversity per sample (sum of branch
  lengths of the minimal subtree spanning each sample's present features). A
  phylogeny-aware alpha-diversity metric that scales to many samples.
- ``tree_layout``: rectangular phylogram coordinates for in-app rendering of
  small trees (capped by tip count; large trees stay download-only).

Pure and dependency-light (Bio.Phylo + the shared matrix loader).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from analysis.sample_comparison import _load_matrix


def _read_tree(tree_path):
    from Bio import Phylo

    return Phylo.read(str(tree_path), "newick")


def _terminal_sets(tree) -> dict:
    """Map each clade -> set of terminal names beneath it (one postorder pass)."""
    sets: dict = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            sets[clade] = {clade.name}
        else:
            union: set = set()
            for child in clade.clades:
                union |= sets.get(child, set())
            sets[clade] = union
    return sets


def _strip_size(feature_id: str) -> str:
    """Drop a trailing VSEARCH ``;size=`` annotation from a feature id."""
    return (feature_id or "").split(";")[0].strip()


def load_asv_mapping(mapping_csv) -> dict:
    """Load ``reversed_asvs.csv`` → ``{original_otu_id (no ;size=): ASV_ID}``.

    The phylo stage renames sequences to ``ASV1..N`` before alignment/tree
    building and records the original headers here, so tree tips can be mapped
    back onto the abundance table's feature ids.
    """
    import csv

    mapping_csv = Path(mapping_csv) if mapping_csv else None
    if not mapping_csv or not mapping_csv.exists() or mapping_csv.stat().st_size == 0:
        return {}
    result: dict = {}
    with mapping_csv.open("r", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle):
            asv = (row.get("ASV_ID") or "").strip()
            original = _strip_size(row.get("Original_Header") or "")
            if asv and original:
                result[original] = asv
    return result


def faith_pd(tree_path, otu_matrix_path, id_mapping_csv=None) -> dict:
    """Per-sample rooted Faith's PD from a Newick tree + feature x sample matrix.

    ``id_mapping_csv`` is the phylo stage's ``reversed_asvs.csv`` mapping tree
    tip ids (``ASV1..N``) to the original OTU/ASV ids used by the abundance
    table. Without it, tips like ``ASV1`` never match feature ids like
    ``A01426:…`` and every sample's PD collapses to zero.
    """
    tree_path = Path(tree_path)
    if not tree_path.exists() or tree_path.stat().st_size == 0:
        return {"error": "Tree not found", "samples": [], "pd": []}
    try:
        tree = _read_tree(tree_path)
    except Exception as exc:  # pragma: no cover - defensive
        return {"error": f"Could not read tree: {exc}", "samples": [], "pd": []}

    feature_ids, samples, values = _load_matrix(Path(otu_matrix_path))
    if not samples:
        return {"error": "No abundance matrix", "samples": [], "pd": []}

    term_sets = _terminal_sets(tree)
    tip_names = term_sets.get(tree.root, set())
    # Precompute (clade, branch_length, terminal_set) for non-root clades once.
    edges = [
        (clade.branch_length or 0.0, term_sets[clade])
        for clade in tree.find_clades()
        if clade is not tree.root
    ]

    # Translate each abundance feature id into the tree's tip namespace. Try, in
    # order: the ASV mapping, the raw id, and the size-stripped id — so this works
    # whether the tree tips are ASV ids or the original headers.
    orig_to_tip = load_asv_mapping(id_mapping_csv)

    def _tip_for(feature_id: str) -> Optional[str]:
        for candidate in (orig_to_tip.get(_strip_size(feature_id)), feature_id, _strip_size(feature_id)):
            if candidate and candidate in tip_names:
                return candidate
        return None

    feature_tip = [_tip_for(fid) for fid in feature_ids]

    matched = 0
    pd_values = []
    for j, sample in enumerate(samples):
        present = {feature_tip[i] for i in range(len(feature_ids)) if values[i][j] > 0 and feature_tip[i]}
        matched = max(matched, len(present))
        total = sum(length for length, tset in edges if tset & present)
        pd_values.append(round(total, 6))

    return {
        "samples": samples,
        "pd": pd_values,
        "total_tips": len(tip_names),
        "matched_features": matched,
    }


def tree_layout(tree_path, max_tips: int = 250) -> Optional[dict]:
    """Rectangular phylogram coordinates for plotting, or ``None`` if too large.

    Returns ``{tips, hlines, vlines}`` where ``tips`` are ``(x, y, name)``,
    ``hlines`` are ``(x_parent, x_child, y)`` and ``vlines`` are ``(x, y0, y1)``.
    """
    tree_path = Path(tree_path)
    if not tree_path.exists() or tree_path.stat().st_size == 0:
        return None
    try:
        tree = _read_tree(tree_path)
    except Exception:  # pragma: no cover - defensive
        return None

    terminals = tree.get_terminals()
    if not terminals or len(terminals) > max_tips:
        return None

    ypos: dict = {}
    for i, tip in enumerate(terminals):
        ypos[tip] = float(i)

    depths = tree.depths()
    if not depths or max(depths.values(), default=0) == 0:
        depths = tree.depths(unit_branch_lengths=True)

    def _y(clade):
        if clade.is_terminal():
            return ypos[clade]
        ys = [_y(c) for c in clade.clades]
        y = sum(ys) / len(ys) if ys else 0.0
        ypos[clade] = y
        return y

    _y(tree.root)

    hlines = []
    vlines = []
    nodes = []
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        x = depths.get(clade, 0.0)
        child_ys = [ypos[c] for c in clade.clades]
        if child_ys:
            vlines.append((x, min(child_ys), max(child_ys)))
        for c in clade.clades:
            hlines.append((x, depths.get(c, 0.0), ypos[c]))
        # Internal-node support, so a viewer can mark weakly supported clades.
        # FastTree parks its local support in ``confidence``; some writers leave
        # it in ``name`` instead.
        support = clade.confidence
        if support is None and clade.name:
            try:
                support = float(clade.name)
            except ValueError:
                support = None
        if support is not None:
            nodes.append((x, ypos[clade], float(support)))

    tips = [(depths.get(t, 0.0), ypos[t], t.name or "") for t in terminals]
    return {"tips": tips, "hlines": hlines, "vlines": vlines, "nodes": nodes, "n_tips": len(terminals)}
