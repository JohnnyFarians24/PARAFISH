"""Persisted taxonomy-stage metrics.

Computed once inside the Databases stage; the Dash callbacks only read the CSVs.

Three things here are easy to get subtly wrong, so they are stated explicitly:

* **Resolution is counted twice** — by feature and by read. Half the OTUs failing
  to reach species means little if they carry 2% of the reads.
* **Contaminants are called from prevalence, not abundance.** A reagent
  contaminant is defined by turning up in the negative controls, not by being
  rare. Samples marked ``is_control`` are never pooled into diversity.
* **Prevalence carries a confidence interval.** "3 of 5 fish infected" is 60%
  with a 95% interval spanning 23%–88%; reporting the point estimate alone
  overstates what five fish can tell you.
"""

from __future__ import annotations

import csv
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from analysis.analysis_metrics import load_matrix

TAXONOMY_RANKS = ["kingdom", "phylum", "class_", "order", "family", "genus", "species"]

_CONSENSUS_COLUMN = {
    "kingdom": "Consensus_Kingdom",
    "phylum": "Consensus_Phylum",
    "class_": "Consensus_Class",
    "order": "Consensus_Order",
    "family": "Consensus_Family",
    "genus": "Consensus_Genus",
    "species": "Consensus_Species",
}

# decontam's prevalence method flags a feature when its presence is significantly
# skewed toward the negative controls. 0.1 is the published default; 0.5 is the
# aggressive setting the decontam vignette recommends when controls are plentiful.
DEFAULT_CONTAMINANT_THRESHOLD = 0.1


def _blank(value) -> bool:
    text = str(value or "").strip().lower()
    return text in {"", "na", "nan", "none", "unclassified", "unknown"}


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _read_consensus(consensus_csv) -> dict[str, dict[str, str]]:
    path = Path(consensus_csv)
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return {row["ASV_ID"]: row for row in csv.DictReader(handle) if row.get("ASV_ID")}


def deepest_consensus_rank(row: dict) -> str:
    """Deepest rank the consensus reached for one feature."""
    for rank in reversed(TAXONOMY_RANKS):
        if not _blank(row.get(_CONSENSUS_COLUMN[rank])):
            return rank
    return "unassigned"


# ── Taxonomic resolution ─────────────────────────────────────────────────


def taxonomic_resolution(consensus_csv, matrix_path=None) -> list[dict]:
    """How many features — and how many reads — reach each taxonomic rank."""
    consensus = _read_consensus(consensus_csv)
    if not consensus:
        return []

    reads_by_feature: dict[str, float] = {}
    if matrix_path and Path(matrix_path).exists():
        features, _, counts = load_matrix(matrix_path)
        totals = counts.sum(axis=0)
        reads_by_feature = {feature: float(totals[index]) for index, feature in enumerate(features)}

    ranks = [*TAXONOMY_RANKS, "unassigned"]
    tally = {rank: {"n_features": 0, "n_reads": 0.0} for rank in ranks}
    for asv_id, row in consensus.items():
        rank = deepest_consensus_rank(row)
        tally[rank]["n_features"] += 1
        tally[rank]["n_reads"] += reads_by_feature.get(asv_id, 0.0)

    total_features = sum(entry["n_features"] for entry in tally.values())
    total_reads = sum(entry["n_reads"] for entry in tally.values())
    return [
        {
            "rank": rank,
            "n_features": tally[rank]["n_features"],
            "pct_features": round(100.0 * tally[rank]["n_features"] / total_features, 3) if total_features else 0.0,
            "n_reads": int(tally[rank]["n_reads"]),
            "pct_reads": round(100.0 * tally[rank]["n_reads"] / total_reads, 3) if total_reads else 0.0,
        }
        for rank in ranks
    ]


# ── Negative controls / decontamination ──────────────────────────────────


def _fisher_greater(present_control: int, absent_control: int, present_sample: int, absent_sample: int) -> float:
    """One-sided Fisher exact p-value for over-representation in the controls."""
    try:
        from scipy.stats import fisher_exact
    except ImportError:  # pragma: no cover - scipy ships with the app
        return 1.0
    _, p_value = fisher_exact([[present_control, absent_control], [present_sample, absent_sample]], alternative="greater")
    return float(p_value)


def decontam_prevalence(
    matrix_path,
    control_samples: list[str],
    threshold: float = DEFAULT_CONTAMINANT_THRESHOLD,
) -> list[dict]:
    """Flag features whose *presence* is skewed toward the negative controls.

    This is decontam's prevalence method: a feature detected in most blanks and
    few real samples is reagent/cross contamination, no matter how abundant it is
    in either. Returns ``[]`` when no controls were marked — guessing which
    samples are blanks would be worse than reporting nothing.
    """
    if not control_samples or not (matrix_path and Path(matrix_path).exists()):
        return []
    features, samples, counts = load_matrix(matrix_path)
    if not features:
        return []

    control_index = [index for index, sample in enumerate(samples) if sample in set(control_samples)]
    true_index = [index for index in range(len(samples)) if index not in control_index]
    if not control_index or not true_index:
        return []

    present = counts > 0
    rows: list[dict] = []
    for column, feature in enumerate(features):
        in_control = int(present[control_index, column].sum())
        in_sample = int(present[true_index, column].sum())
        score = _fisher_greater(in_control, len(control_index) - in_control, in_sample, len(true_index) - in_sample)
        rows.append(
            {
                "feature": feature,
                "prevalence_controls": round(in_control / len(control_index), 4),
                "prevalence_samples": round(in_sample / len(true_index), 4),
                "reads_controls": int(counts[control_index, column].sum()),
                "reads_samples": int(counts[true_index, column].sum()),
                "score": round(score, 6),
                "is_contaminant": int(score < threshold),
            }
        )
    rows.sort(key=lambda row: row["score"])
    return rows


# ── Prevalence with a confidence interval ────────────────────────────────


def wilson_interval(successes: int, trials: int, z: float = 1.96) -> tuple[float, float]:
    """95% Wilson score interval for a proportion.

    Preferred over the normal approximation, which produces intervals outside
    [0, 1] — and negative lower bounds — exactly where parasitology lives: few
    hosts and prevalences near 0 or 1.
    """
    if trials <= 0:
        return 0.0, 0.0
    phat = successes / trials
    denominator = 1 + z**2 / trials
    centre = phat + z**2 / (2 * trials)
    spread = z * math.sqrt(phat * (1 - phat) / trials + z**2 / (4 * trials**2))
    return max(0.0, (centre - spread) / denominator), min(1.0, (centre + spread) / denominator)


def taxon_prevalence(matrix_path, consensus_csv, rank: str = "genus", exclude_samples: Optional[list[str]] = None) -> list[dict]:
    """Per-taxon prevalence across biological samples, with a 95% interval."""
    if not (matrix_path and Path(matrix_path).exists()):
        return []
    consensus = _read_consensus(consensus_csv)
    features, samples, counts = load_matrix(matrix_path)
    if not features:
        return []

    excluded = set(exclude_samples or [])
    keep = [index for index, sample in enumerate(samples) if sample not in excluded]
    if not keep:
        return []
    counts = counts[keep]
    n_samples = len(keep)

    column_name = _CONSENSUS_COLUMN.get(rank, "Consensus_Genus")
    by_taxon: dict[str, np.ndarray] = {}
    for column, feature in enumerate(features):
        label = (consensus.get(feature) or {}).get(column_name, "")
        taxon = "Unassigned" if _blank(label) else str(label).strip()
        vector = counts[:, column]
        by_taxon[taxon] = by_taxon[taxon] + vector if taxon in by_taxon else vector.copy()

    rows = []
    for taxon, vector in by_taxon.items():
        detected = int((vector > 0).sum())
        low, high = wilson_interval(detected, n_samples)
        rows.append(
            {
                "taxon": taxon,
                "n_detected": detected,
                "n_samples": n_samples,
                "prevalence": round(detected / n_samples, 4),
                "ci_low": round(low, 4),
                "ci_high": round(high, 4),
                "total_reads": int(vector.sum()),
            }
        )
    rows.sort(key=lambda row: (-row["prevalence"], -row["total_reads"]))
    return rows


# ── Grouped beta diversity ───────────────────────────────────────────────


def permanova_by_group(matrix_path, groups: dict[str, str], permutations: int = 999) -> dict:
    """PERMANOVA on Bray-Curtis distances between the labelled sample groups.

    Returns ``{}`` when there are not at least two groups with two samples each —
    PERMANOVA on a singleton group tests nothing.
    """
    if not (matrix_path and Path(matrix_path).exists()) or not groups:
        return {}
    _, samples, counts = load_matrix(matrix_path)
    keep = [index for index, sample in enumerate(samples) if groups.get(sample)]
    if len(keep) < 3:
        return {}
    labels = [groups[samples[index]] for index in keep]
    sizes: dict[str, int] = {}
    for label in labels:
        sizes[label] = sizes.get(label, 0) + 1
    if len(sizes) < 2 or min(sizes.values()) < 2:
        return {}

    from skbio.stats.distance import permanova
    from skbio.diversity import beta_diversity

    ids = [samples[index] for index in keep]
    distances = beta_diversity("braycurtis", counts[keep].astype(int), ids)
    result = permanova(distances, labels, permutations=permutations)
    return {
        "test": "PERMANOVA (Bray-Curtis)",
        "n_samples": len(ids),
        "n_groups": len(sizes),
        "pseudo_f": round(float(result["test statistic"]), 4),
        "p_value": float(result["p-value"]),
        "permutations": permutations,
    }


# ── Exports: taxa tables at every rank, and BIOM ─────────────────────────


def taxa_tables_all_ranks(matrix_path, consensus_csv, out_dir) -> dict[str, Path]:
    """Write a taxon x sample table for *every* rank, absolute and relative.

    Only a genus table used to be produced, which forced anyone wanting phylum- or
    species-level composition to redo the collapse by hand.
    """
    from analysis.community_analysis import collapsed_taxa_table

    if not (matrix_path and Path(matrix_path).exists() and consensus_csv and Path(consensus_csv).exists()):
        return {}
    out_dir = Path(out_dir)
    written: dict[str, Path] = {}
    for rank in TAXONOMY_RANKS:
        absolute = collapsed_taxa_table(matrix_path, consensus_csv, rank, out_dir / f"taxa_table_{rank.rstrip('_')}.csv")
        if not absolute:
            continue
        written[f"taxa_table_{rank.rstrip('_')}_csv"] = Path(absolute)
        relative = _relative_taxa_table(Path(absolute), out_dir / f"taxa_table_{rank.rstrip('_')}_relative.csv")
        if relative:
            written[f"taxa_table_{rank.rstrip('_')}_relative_csv"] = relative
    return written


def _relative_taxa_table(absolute_csv: Path, out_path: Path) -> Optional[Path]:
    """Convert a taxon x sample count table to per-sample percentages."""
    with absolute_csv.open("r", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 2:
        return None
    header, *body = rows
    totals = [0.0] * (len(header) - 1)
    for row in body:
        for index, cell in enumerate(row[1:]):
            totals[index] += float(cell or 0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in body:
            values = [
                round(100.0 * float(cell or 0) / totals[index], 6) if totals[index] else 0.0
                for index, cell in enumerate(row[1:])
            ]
            writer.writerow([row[0], *values])
    return out_path


def write_biom(matrix_path, consensus_csv, out_path) -> Optional[Path]:
    """Write a BIOM 1.0 (JSON) table with taxonomy metadata per feature.

    Hand-rolled rather than pulled in via the ``biom-format`` package: the sparse
    JSON schema is small, and this keeps the export dependency-free and readable.
    """
    import json

    if not (matrix_path and Path(matrix_path).exists()):
        return None
    features, samples, counts = load_matrix(matrix_path)
    if not features or not samples:
        return None
    consensus = _read_consensus(consensus_csv) if consensus_csv else {}

    def _lineage(feature: str) -> list[str]:
        row = consensus.get(feature) or {}
        prefixes = {"kingdom": "k__", "phylum": "p__", "class_": "c__", "order": "o__", "family": "f__", "genus": "g__", "species": "s__"}
        lineage = []
        for rank in TAXONOMY_RANKS:
            value = row.get(_CONSENSUS_COLUMN[rank])
            lineage.append(f"{prefixes[rank]}{'' if _blank(value) else str(value).strip()}")
        return lineage

    # BIOM's sparse "data" is [row, column, value] over non-zero cells only.
    data = [
        [feature_index, sample_index, int(counts[sample_index, feature_index])]
        for feature_index in range(len(features))
        for sample_index in range(len(samples))
        if counts[sample_index, feature_index] > 0
    ]
    table = {
        "id": Path(matrix_path).stem,
        "format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "PARAFISH",
        "date": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [len(features), len(samples)],
        "rows": [{"id": feature, "metadata": {"taxonomy": _lineage(feature)}} for feature in features],
        "columns": [{"id": sample, "metadata": None} for sample in samples],
        "data": data,
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(table), encoding="utf-8")
    return out_path


_PHYLOSEQ_README = '''# Load this PARAFISH run into a phyloseq object.
#
# Run this from the directory that holds it — the three tables it reads are its
# neighbours. Written by PARAFISH; regenerated on every Databases run.
#
#   otu_table.tsv    feature x sample counts (raw, not rarefied, not relative)
#   tax_table.tsv    feature -> Kingdom..Species, blank where the consensus stopped
#   sample_data.tsv  sample -> group / control flag / read total
#
# Why this exists next to table.biom: BIOM carries the taxonomy but PARAFISH's BIOM
# has no per-sample metadata, so the grouping and the control flags — the two things
# every comparison needs — do not survive that route. They do survive this one.

library(phyloseq)

otu <- as.matrix(read.delim("otu_table.tsv", row.names = 1, check.names = FALSE))
tax <- as.matrix(read.delim("tax_table.tsv", row.names = 1, na.strings = ""))
sam <- read.delim("sample_data.tsv", row.names = 1, check.names = FALSE)

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  tax_table(tax),
  sample_data(sam)
)

# Negative/blank controls are kept in the object rather than dropped, so that the
# decision to exclude them stays yours and stays visible:
#   ps_samples <- prune_samples(sample_data(ps)$is_control == "no", ps)

# The tree in phylo/<feature_set>/tree.nwk is NOT merged in here: its tips were
# renamed to ASV1..N during alignment. Map them back with the id table first —
# phylo/mafft/reversed_asvs.csv — or the tips will not match these feature ids.

ps
'''


def write_phyloseq_bundle(
    matrix_path,
    consensus_csv,
    out_dir,
    *,
    control_samples: Optional[list[str]] = None,
    groups: Optional[dict[str, str]] = None,
) -> dict[str, Path]:
    """Write the three tables phyloseq imports, plus the R script that assembles them.

    This is not redundant with :func:`write_biom`. BIOM's sample column carries
    ``metadata: None`` here, so a run exported that way arrives in R with its groups
    and its control flags gone — and those are exactly what the downstream comparison
    is about. Splitting the export into phyloseq's own three tables keeps them.

    Counts are raw: rarefying or converting to relative abundance is a decision that
    belongs to the analysis, and phyloseq has ``rarefy_even_depth``/``transform_sample_counts``
    for it. Exporting pre-transformed counts would quietly make that choice for the user.
    """
    if not (matrix_path and Path(matrix_path).exists()):
        return {}
    features, samples, counts = load_matrix(matrix_path)
    if not features or not samples:
        return {}

    consensus = _read_consensus(consensus_csv) if consensus_csv else {}
    controls = {str(sample) for sample in (control_samples or [])}
    groups = groups or {}
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    otu_path = out_dir / "otu_table.tsv"
    with otu_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["feature_id", *samples])
        for feature_index, feature in enumerate(features):
            writer.writerow([feature, *(int(counts[sample_index, feature_index]) for sample_index in range(len(samples)))])
    written["phyloseq_otu_table"] = otu_path

    # phyloseq wants one column per rank, and a *blank* where the assignment stopped.
    # Writing "unclassified" instead would invent a taxon that tax_glom() then treats
    # as a real one, silently merging unrelated features under it.
    tax_path = out_dir / "tax_table.tsv"
    with tax_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["feature_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])
        for feature in features:
            row = consensus.get(feature) or {}
            values = []
            for rank in TAXONOMY_RANKS:
                value = row.get(_CONSENSUS_COLUMN[rank])
                values.append("" if _blank(value) else str(value).strip())
            writer.writerow([feature, *values])
    written["phyloseq_tax_table"] = tax_path

    sample_path = out_dir / "sample_data.tsv"
    with sample_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "group", "is_control", "reads"])
        for sample_index, sample in enumerate(samples):
            writer.writerow(
                [
                    sample,
                    groups.get(sample, "") or "ungrouped",
                    "yes" if sample in controls else "no",
                    int(counts[sample_index].sum()),
                ]
            )
    written["phyloseq_sample_data"] = sample_path

    script_path = out_dir / "load_phyloseq.R"
    script_path.write_text(_PHYLOSEQ_README, encoding="utf-8")
    written["phyloseq_script"] = script_path
    return written


# ── Persistence ──────────────────────────────────────────────────────────


def build_database_metrics(
    databases_dir,
    consensus_csv,
    matrix_path=None,
    *,
    control_samples: Optional[list[str]] = None,
    groups: Optional[dict[str, str]] = None,
    rank: str = "genus",
) -> dict[str, Path]:
    """Compute and persist the taxonomy-stage tables. Returns name -> path."""
    databases_dir = Path(databases_dir)
    metrics_dir = databases_dir / "metrics"
    written: dict[str, Path] = {}

    if not (consensus_csv and Path(consensus_csv).exists()):
        return written

    resolution = taxonomic_resolution(consensus_csv, matrix_path)
    if resolution:
        written["resolution_csv"] = _write_csv(
            metrics_dir / "taxonomic_resolution.csv",
            ["rank", "n_features", "pct_features", "n_reads", "pct_reads"],
            resolution,
        )

    control_samples = control_samples or []
    contaminants = decontam_prevalence(matrix_path, control_samples)
    if contaminants:
        written["contaminants_csv"] = _write_csv(
            metrics_dir / "contaminants.csv",
            ["feature", "prevalence_controls", "prevalence_samples", "reads_controls", "reads_samples", "score", "is_contaminant"],
            contaminants,
        )

    prevalence = taxon_prevalence(matrix_path, consensus_csv, rank=rank, exclude_samples=control_samples)
    if prevalence:
        written["taxon_prevalence_csv"] = _write_csv(
            metrics_dir / "taxon_prevalence.csv",
            ["taxon", "n_detected", "n_samples", "prevalence", "ci_low", "ci_high", "total_reads"],
            prevalence,
        )

    if groups:
        try:
            result = permanova_by_group(matrix_path, groups)
        except Exception:
            result = {}
        if result:
            written["permanova_csv"] = _write_csv(metrics_dir / "permanova.csv", list(result), [result])

    # Taxon x sample tables at every rank (absolute and relative), plus a BIOM
    # export so the run can be taken straight into phyloseq/QIIME.
    try:
        written.update(taxa_tables_all_ranks(matrix_path, consensus_csv, databases_dir / "taxa_tables"))
    except Exception:
        pass
    try:
        biom = write_biom(matrix_path, consensus_csv, databases_dir / "table.biom")
        if biom:
            written["biom_json"] = biom
    except Exception:
        pass
    # BIOM loses the sample groups and control flags on the way into R; the phyloseq
    # bundle is the route that keeps them.
    try:
        written.update(
            write_phyloseq_bundle(
                matrix_path,
                consensus_csv,
                databases_dir / "phyloseq",
                control_samples=control_samples,
                groups=groups,
            )
        )
    except Exception:
        pass

    return written


def database_metrics_paths(databases_dir, feature_set: str = "otu") -> dict[str, Path]:
    """Return the persisted taxonomy tables that exist under the databases STAGE dir.

    ``databases_dir`` is ``<run>/databases``; the tables live one level down, under the
    feature set the run produced (``databases/otu/metrics/``). This mirrors
    ``analysis.analysis_metrics.analysis_metrics_paths`` deliberately — the two stages
    have the same layout and had two different opinions about it, which is how the
    Results page came to read a directory the pipeline never wrote.
    """
    names = {
        "resolution_csv": "taxonomic_resolution.csv",
        "contaminants_csv": "contaminants.csv",
        "taxon_prevalence_csv": "taxon_prevalence.csv",
        "permanova_csv": "permanova.csv",
    }
    # Runs predating feature sets wrote straight into databases/metrics/.
    for metrics_dir in (Path(databases_dir) / feature_set / "metrics", Path(databases_dir) / "metrics"):
        found = {key: metrics_dir / name for key, name in names.items() if (metrics_dir / name).exists()}
        if found:
            return found
    return {}


# The taxonomy stage's per-feature-set tables, by filename. These sit beside the
# metrics/ dir rather than inside it.
_TAXONOMY_TABLES = {
    "taxonomy_consensus.csv",
    "source_assignments.csv",
    "source_status.csv",
    "conflict_report.csv",
    "taxonomy_hits.csv",
}


def taxonomy_table_path(databases_dir, filename: str, feature_set: str = "otu") -> Optional[Path]:
    """Path to one of the taxonomy stage's tables, or ``None``.

    ``ui.feature_sets.taxonomy_file`` does this for the UI layer, resolving the store
    first; this is the plain on-disk half, so non-UI readers (``analysis.results_summary``)
    can resolve the same layout without importing from ``ui``.
    """
    for candidate in (Path(databases_dir) / feature_set / filename, Path(databases_dir) / filename):
        if candidate.exists():
            return candidate
    return None
