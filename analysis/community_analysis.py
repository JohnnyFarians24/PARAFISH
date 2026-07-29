"""Community ecology outputs for cross-sample comparison at scale.

Turns the OTU/ASV abundance matrix (+ the taxonomy consensus table) into the
outputs metabarcoding researchers actually use:

- taxonomic composition (relative abundance, collapsed to a rank, top-N + Other);
- a collapsed taxon x sample table;
- per-sample alpha diversity + a ranking summary (find the best/worst samples);
- beta diversity (Bray-Curtis) and a PCoA ordination so 100 samples reduce to a
  2D map where similar communities cluster.

Pure and numpy-only (no scipy/pandas) so it stays unit-testable and the UI layer
turns the returned dicts into figures.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Optional

import numpy as np

from analysis.sample_comparison import _load_matrix


# Consensus column candidates per rank (taxonomy_consensus.csv uses Consensus_*).
_RANK_COLUMNS = {
    "kingdom": ["Consensus_Kingdom", "kingdom", "Kingdom"],
    "phylum": ["Consensus_Phylum", "phylum", "Phylum"],
    "class": ["Consensus_Class", "class", "class_", "Class"],
    "order": ["Consensus_Order", "order", "Order"],
    "family": ["Consensus_Family", "family", "Family"],
    "genus": ["Consensus_Genus", "genus", "Genus"],
    "species": ["Consensus_Species", "species", "Species"],
}

# ``database_metrics.TAXONOMY_RANKS`` spells this rank "class_" (dodging the Python
# keyword) and passes that string straight into the functions below, while every other
# rank has a single spelling. The lookup is a plain dict get, so the mismatch never
# raised: it returned {}, and taxa_table_class.csv reported 100% "Unclassified" for a
# run whose tax_table said Craniata. One rank, one silent hole, because one module
# spells it with an underscore. Accepted as an alias rather than duplicated into
# _RANK_COLUMNS, so RANKS below stays the canonical seven.
_RANK_ALIASES = {"class_": "class"}


def canonical_rank(rank: str) -> str:
    """The one spelling of ``rank`` this module keys on."""
    rank = str(rank or "").strip()
    return _RANK_ALIASES.get(rank, rank)

RANKS = list(_RANK_COLUMNS.keys())

_UNCLASSIFIED = {"", "na", "nan", "none", "unclassified", "unassigned", "unknown"}


def _norm_taxon(value: Optional[str]) -> str:
    text = (value or "").strip()
    return "" if text.lower() in _UNCLASSIFIED else text


def _load_taxonomy(path: Path, rank: str) -> dict[str, str]:
    """Map feature id -> taxon name at ``rank`` from a consensus CSV."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    rows = list(csv.DictReader(text.splitlines()))
    if not rows:
        return {}
    columns = list(rows[0].keys())
    id_col = next((c for c in ("feature_id", "ASV_ID", "asv_id", "OTU_ID", "id", "qseqid") if c in columns), columns[0])
    rank_col = next((c for c in _RANK_COLUMNS.get(canonical_rank(rank), []) if c in columns), None)
    if rank_col is None:
        return {}
    mapping: dict[str, str] = {}
    for row in rows:
        feature = (row.get(id_col) or "").strip()
        if feature:
            mapping[feature] = _norm_taxon(row.get(rank_col))
    return mapping


def taxa_composition(otu_matrix_path, taxonomy_csv_path, rank: str = "genus", top_n: int = 12) -> dict:
    """Collapse the abundance matrix by taxon at ``rank`` and return relative
    composition with the ``top_n`` most abundant taxa (rest folded into 'Other').
    """
    feature_ids, samples, values = _load_matrix(Path(otu_matrix_path))
    if not samples or not feature_ids:
        return {"error": "No abundance matrix", "samples": [], "taxa": [], "rel": [], "abs": []}
    taxonomy = _load_taxonomy(Path(taxonomy_csv_path), rank) if taxonomy_csv_path else {}

    # Sum abundances per taxon across features.
    taxon_totals: dict[str, list[float]] = {}
    for i, feature in enumerate(feature_ids):
        taxon = taxonomy.get(feature, "") or "Unclassified"
        if taxon not in taxon_totals:
            taxon_totals[taxon] = [0.0] * len(samples)
        row = values[i]
        acc = taxon_totals[taxon]
        for j in range(len(samples)):
            acc[j] += row[j]

    # Rank taxa by grand total; keep top_n, fold the rest into "Other".
    ranked = sorted(taxon_totals.items(), key=lambda kv: sum(kv[1]), reverse=True)
    keep = ranked[: max(1, top_n)]
    rest = ranked[max(1, top_n):]
    taxa = [name for name, _ in keep]
    abs_rows = [counts[:] for _, counts in keep]
    if rest:
        other = [0.0] * len(samples)
        for _, counts in rest:
            for j in range(len(samples)):
                other[j] += counts[j]
        taxa.append("Other")
        abs_rows.append(other)

    # Relative abundance per sample (column-normalised).
    col_totals = [sum(abs_rows[t][j] for t in range(len(taxa))) for j in range(len(samples))]
    rel_rows = [
        [(abs_rows[t][j] / col_totals[j] if col_totals[j] > 0 else 0.0) for j in range(len(samples))]
        for t in range(len(taxa))
    ]
    return {"rank": rank, "samples": samples, "taxa": taxa, "rel": rel_rows, "abs": abs_rows}


def collapsed_taxa_table(otu_matrix_path, taxonomy_csv_path, rank: str, out_path) -> Optional[Path]:
    """Write a taxon x sample absolute-abundance CSV; return the path."""
    comp = taxa_composition(otu_matrix_path, taxonomy_csv_path, rank=rank, top_n=10**9)
    if comp.get("error"):
        return None
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([rank] + comp["samples"])
        for t, taxon in enumerate(comp["taxa"]):
            writer.writerow([taxon] + [int(round(v)) for v in comp["abs"][t]])
    return out_path


def _shannon(counts: list[float]) -> float:
    total = sum(c for c in counts if c > 0)
    if total <= 0:
        return 0.0
    return -sum((c / total) * math.log(c / total) for c in counts if c > 0)


def sample_summary(otu_matrix_path, taxonomy_csv_path=None, rank: str = "genus") -> list[dict]:
    """Per-sample summary for ranking: reads, richness, Shannon, % classified."""
    feature_ids, samples, values = _load_matrix(Path(otu_matrix_path))
    if not samples:
        return []
    taxonomy = _load_taxonomy(Path(taxonomy_csv_path), rank) if taxonomy_csv_path else {}
    summary = []
    for j, sample in enumerate(samples):
        counts = [values[i][j] for i in range(len(feature_ids))]
        present = [i for i in range(len(feature_ids)) if counts[i] > 0]
        classified_reads = sum(counts[i] for i in present if taxonomy.get(feature_ids[i]))
        total_reads = sum(counts)
        summary.append(
            {
                "sample": sample,
                "reads": int(total_reads),
                "richness": len(present),
                "shannon": round(_shannon(counts), 4),
                "pct_classified": round(100.0 * classified_reads / total_reads, 1) if total_reads > 0 else 0.0,
            }
        )
    return summary


def bray_curtis_matrix(otu_matrix_path) -> dict:
    """Sample x sample Bray-Curtis dissimilarity (0 = identical, 1 = disjoint)."""
    feature_ids, samples, values = _load_matrix(Path(otu_matrix_path))
    if len(samples) < 2:
        return {"samples": samples, "matrix": [[0.0]] if samples else []}
    mat = np.array(values, dtype=float)  # features x samples
    n = len(samples)
    bc = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            x, y = mat[:, a], mat[:, b]
            denom = (x + y).sum()
            d = float(np.abs(x - y).sum() / denom) if denom > 0 else 0.0
            bc[a, b] = bc[b, a] = d
    return {"samples": samples, "matrix": bc.tolist()}


def pcoa(dissimilarity: dict) -> dict:
    """Principal coordinates analysis (classical MDS) of a dissimilarity matrix.

    Returns per-sample 2D coordinates and the variance explained by PC1/PC2.
    Uses numpy only (double-centering + symmetric eigendecomposition).
    """
    samples = dissimilarity.get("samples") or []
    matrix = dissimilarity.get("matrix") or []
    n = len(samples)
    if n < 3:
        # Not enough points for a meaningful 2D ordination.
        return {"samples": samples, "coords": [[0.0, 0.0] for _ in samples], "variance": [0.0, 0.0]}
    d = np.array(matrix, dtype=float)
    a = -0.5 * (d ** 2)
    j = np.eye(n) - np.ones((n, n)) / n
    b = j @ a @ j
    eigvals, eigvecs = np.linalg.eigh(b)  # ascending
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals[eigvals > 0].sum()
    coords = []
    for i in range(n):
        pc1 = eigvecs[i, 0] * math.sqrt(eigvals[0]) if eigvals[0] > 0 else 0.0
        pc2 = eigvecs[i, 1] * math.sqrt(eigvals[1]) if len(eigvals) > 1 and eigvals[1] > 0 else 0.0
        coords.append([float(pc1), float(pc2)])
    variance = [
        round(100.0 * eigvals[0] / pos, 1) if pos > 0 and eigvals[0] > 0 else 0.0,
        round(100.0 * eigvals[1] / pos, 1) if pos > 0 and len(eigvals) > 1 and eigvals[1] > 0 else 0.0,
    ]
    return {"samples": samples, "coords": coords, "variance": variance}
