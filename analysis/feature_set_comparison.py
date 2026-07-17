"""Compare the OTU and ASV feature sets without ever merging them.

An OTU is a cluster of reads within 97 % identity of a centroid; an ASV is an exact
denoised sequence. One 97 % OTU routinely contains several ASVs, so their counts are
not on the same scale and their ids mean different things. Concatenating the two
tables — or joining them on feature id — produces a table that is wrong in a way
that is very hard to see afterwards.

The two levels at which they CAN legitimately be compared are:

* **Summary** — the same statistic computed on each set, reported side by side, one
  row per set. Never one row containing both.
* **Taxonomy** — both sets classify against the same references, so "did OTUs and
  ASVs both detect *Diplostomum*, and with how many features and reads each?" is a
  real question with a real answer. The bridge is the taxon name, not the feature id.

Everything here writes SEPARATE tables and joins only on taxon.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Optional

from analysis.analysis_metrics import load_matrix

SUMMARY_COLUMNS = [
    "feature_set",
    "label",
    "n_features",
    "n_samples",
    "total_reads",
    "features_in_multiple_samples",
    "features_in_one_sample",
    "mean_reads_per_feature",
]

TAXON_COLUMNS = [
    "taxon",
    "rank",
    "in_otu",
    "in_asv",
    "n_otus",
    "n_asvs",
    "otu_reads",
    "asv_reads",
    "agreement",
]

_RANK_COLUMN = {
    "species": "Consensus_Species",
    "genus": "Consensus_Genus",
    "family": "Consensus_Family",
}


def _read_consensus(path: Optional[Path]) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle))


def _reads_by_feature(matrix_path: Optional[Path]) -> dict[str, float]:
    if not matrix_path or not Path(matrix_path).exists():
        return {}
    features, _samples, counts = load_matrix(matrix_path)
    if not features:
        return {}
    return {feature: float(counts[:, index].sum()) for index, feature in enumerate(features)}


def build_summary_rows(feature_sets: list[dict]) -> list[dict]:
    """One row per feature set — the same numbers, computed the same way, side by side."""
    rows = []
    for feature_set in feature_sets:
        matrix = feature_set.get("matrix")
        n_features = n_samples = 0
        total_reads = 0.0
        shared = unique = 0
        if matrix and Path(matrix).exists():
            features, samples, counts = load_matrix(matrix)
            n_features = len(features)
            n_samples = len(samples)
            total_reads = float(counts.sum())
            present = counts > 0
            for index in range(n_features):
                hits = int(present[:, index].sum())
                if hits > 1:
                    shared += 1
                elif hits == 1:
                    unique += 1
        rows.append(
            {
                "feature_set": feature_set["key"],
                "label": feature_set["label"],
                "n_features": n_features,
                "n_samples": n_samples,
                "total_reads": int(total_reads),
                "features_in_multiple_samples": shared,
                "features_in_one_sample": unique,
                "mean_reads_per_feature": round(total_reads / n_features, 2) if n_features else 0.0,
            }
        )
    return rows


def build_taxon_bridge(
    otu_consensus: Optional[Path],
    asv_consensus: Optional[Path],
    otu_matrix: Optional[Path],
    asv_matrix: Optional[Path],
    rank: str = "genus",
) -> list[dict]:
    """Taxa found by each feature set, and whether the two agree.

    This is the only legitimate join between the sets: on the taxon name. Each side
    keeps its own feature count and its own read count — they are reported in
    separate columns and never summed together.
    """
    column = _RANK_COLUMN.get(rank, "Consensus_Genus")

    def _tally(consensus_path: Optional[Path], matrix_path: Optional[Path]) -> dict[str, dict]:
        reads = _reads_by_feature(matrix_path)
        tally: dict[str, dict] = defaultdict(lambda: {"n": 0, "reads": 0.0})
        for row in _read_consensus(consensus_path):
            taxon = (row.get(column) or "").strip()
            if not taxon:
                continue
            feature_id = (row.get("ASV_ID") or "").strip()
            tally[taxon]["n"] += 1
            tally[taxon]["reads"] += reads.get(feature_id, 0.0)
        return tally

    otu = _tally(otu_consensus, otu_matrix)
    asv = _tally(asv_consensus, asv_matrix)

    rows = []
    for taxon in sorted(set(otu) | set(asv)):
        in_otu = taxon in otu
        in_asv = taxon in asv
        rows.append(
            {
                "taxon": taxon,
                "rank": rank,
                "in_otu": int(in_otu),
                "in_asv": int(in_asv),
                "n_otus": otu.get(taxon, {}).get("n", 0),
                "n_asvs": asv.get(taxon, {}).get("n", 0),
                "otu_reads": int(otu.get(taxon, {}).get("reads", 0)),
                "asv_reads": int(asv.get(taxon, {}).get("reads", 0)),
                "agreement": "both" if in_otu and in_asv else ("otu_only" if in_otu else "asv_only"),
            }
        )
    # Taxa both methods agree on first, then the biggest disagreements: a taxon only
    # one method sees is exactly what a reviewer will ask about.
    rows.sort(key=lambda row: (row["agreement"] != "both", -(row["otu_reads"] + row["asv_reads"])))
    return rows


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def build_feature_set_comparison(
    output_dir,
    feature_sets: list[dict],
    consensus_by_set: Optional[dict] = None,
) -> dict[str, Path]:
    """Write the OTU-vs-ASV comparison tables. Returns name -> path."""
    output_dir = Path(output_dir)
    written: dict[str, Path] = {}
    if len(feature_sets) < 2:
        return written

    written["feature_set_summary_csv"] = _write_csv(
        output_dir / "feature_set_summary.csv", SUMMARY_COLUMNS, build_summary_rows(feature_sets)
    )

    by_key = {feature_set["key"]: feature_set for feature_set in feature_sets}
    consensus_by_set = consensus_by_set or {}
    bridge = build_taxon_bridge(
        consensus_by_set.get("otu"),
        consensus_by_set.get("asv"),
        by_key.get("otu", {}).get("matrix"),
        by_key.get("asv", {}).get("matrix"),
    )
    if bridge:
        written["feature_set_taxa_csv"] = _write_csv(output_dir / "feature_set_taxa.csv", TAXON_COLUMNS, bridge)

    return written
