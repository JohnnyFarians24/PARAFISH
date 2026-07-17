"""Persisted analysis metrics: chimeras, rarefaction, alpha diversity, rank-abundance.

Everything here runs once, inside the analysis stage, and writes CSVs the Dash
callbacks read back. Rarefaction in particular is a Monte-Carlo subsampling over a
grid of depths and must never run per render.

Two conventions the rest of the app depends on:

* Diversity is computed on a **rarefied** matrix — every sample subsampled to the
  same depth — because richness grows with sequencing effort. The depth used is
  written into ``rarefaction_depth.csv`` so the UI can state it.
* Chimeras are reported twice: as a share of **variants** and as a share of
  **reads**. A run can have 40% chimeric variants that carry 0.5% of the reads;
  quoting either number alone is misleading.
"""

from __future__ import annotations

import csv
import math
import re
from collections import Counter
from pathlib import Path
from typing import Iterator, Optional

import numpy as np

# Depths sampled along the rarefaction curve, and how many random subsamples are
# averaged at each depth. 10 iterations is enough to smooth the curve without
# making the stage crawl on runs with hundreds of samples.
RAREFACTION_POINTS = 20
RAREFACTION_ITERATIONS = 10

# A feature present in at least this share of samples is "core".
DEFAULT_CORE_FRACTION = 0.8

ALPHA_COLUMNS = ["sample", "reads", "observed", "chao1", "shannon", "simpson", "pielou", "hill_q0", "hill_q1", "hill_q2"]

# The two feature sets the platform can produce. They answer different questions
# and are never pooled into one table: an OTU is a 97%-identity cluster, an ASV is
# an exact denoised sequence, so their counts are not on the same scale.
FEATURE_SET_LABELS = {"otu": "OTU (VSEARCH)", "asv": "ASV (DADA2)"}


# ── FASTA / chimera accounting ───────────────────────────────────────────

_SIZE = re.compile(r"size=(\d+)")


def _fasta_headers(path: Path) -> Iterator[str]:
    with path.open("r", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip()


def fasta_size_total(path) -> tuple[int, int]:
    """Return ``(n_sequences, total_reads)`` for a VSEARCH ``;size=``-annotated FASTA.

    Sequences with no ``size=`` count as one read, which is what VSEARCH assumes.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return 0, 0
    count = 0
    reads = 0
    for header in _fasta_headers(path):
        count += 1
        match = _SIZE.search(header)
        reads += int(match.group(1)) if match else 1
    return count, reads


def fasta_lengths(path) -> Counter:
    """Length histogram of a FASTA, streaming."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return Counter()
    counts: Counter = Counter()
    length = 0
    with path.open("r", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if length:
                    counts[length] += 1
                length = 0
            else:
                length += len(line.strip())
    if length:
        counts[length] += 1
    return counts


def original_label(header: str) -> str:
    """The pre-relabel header VSEARCH keeps after the new id (``--relabel_keep``).

    The chimera step relabels its output to OTU1..n, and — importantly — restarts
    that numbering in the *chimeras* file too, so "OTU1" means different sequences
    in the two outputs. Only the kept original label is a stable key that also
    matches the cluster map, which is written before relabelling.
    """
    header = header.strip()
    parts = header.split(None, 1)
    return parts[1].strip() if len(parts) > 1 else header


def cluster_read_counts(clusters_uc) -> dict[str, int]:
    """Reads behind each cluster centroid, from a VSEARCH ``--uc`` file.

    The ``.uc`` file is written by the clustering step, before relabelling, so its
    keys are the original centroid headers. ``S`` rows name the centroid, ``H``
    rows name each member and the centroid it joined.
    """
    clusters_uc = Path(clusters_uc)
    if not clusters_uc.exists():
        return {}

    def _size(label: str) -> int:
        match = _SIZE.search(label)
        return int(match.group(1)) if match else 1

    totals: dict[str, int] = {}
    with clusters_uc.open("r", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            record, query, target = fields[0], fields[8], fields[9]
            if record == "S":
                totals[query] = totals.get(query, 0) + _size(query)
            elif record == "H":
                totals[target] = totals.get(target, 0) + _size(query)
    return totals


def matrix_read_total(matrix_path) -> Optional[int]:
    """The run's canonical read total: the sum of the mapping-based abundance table."""
    if not (matrix_path and Path(matrix_path).exists()):
        return None
    _, samples, counts = load_matrix(matrix_path)
    if not samples:
        return None
    return int(counts.sum())


def chimera_stats(chimeras_fasta, otus_raw_fasta, clusters_uc=None, matrix_path=None) -> dict:
    """Chimera rate as a share of variants *and* as a share of reads.

    The two numbers answer different questions and neither implies the other: a run
    can flag 40% of the variants as chimeric while they carry 0.5% of the reads.

    ``total_reads`` is the canonical abundance total — the sum of the mapping-based
    ``otu_table``, the same number every other table in the run reports. It used to be
    ``sum(cluster_read_counts())`` instead, which is a *different population*: the
    clustering never sees the reads that ``--minuniquesize`` dropped at dereplication
    (10 527 of them on the reference run), so this file said 147 010 while otu_table,
    taxonomic_resolution and every taxa_table said 157 527.

    ``chimeric_reads`` necessarily still comes from the cluster map, and that is not an
    inconsistency: chimeric variants are removed *before* the mapping step, so their
    reads are absent from otu_table by construction — the cluster map is the only place
    they are counted at all. The percentage is therefore taken over
    ``chimeric + mapped``: every read that had a chance to be counted.
    """
    n_chimeric, chimeric_size = fasta_size_total(chimeras_fasta)
    n_total, total_size = fasta_size_total(otus_raw_fasta)

    cluster_reads = cluster_read_counts(clusters_uc) if clusters_uc else {}
    if cluster_reads:
        # Match on the kept original label: after --relabel the chimera file's own
        # "OTU1" is a different sequence from the non-chimera file's "OTU1", and
        # neither appears in the cluster map.
        chimeric_labels = (
            {original_label(header) for header in _fasta_headers(Path(chimeras_fasta))}
            if Path(chimeras_fasta).exists()
            else set()
        )
        chimeric_reads = sum(reads for label, reads in cluster_reads.items() if label in chimeric_labels)
    else:
        chimeric_reads = chimeric_size

    mapped_reads = matrix_read_total(matrix_path)
    if mapped_reads is None:
        # No abundance table (a failed or partial run): fall back to what the clustering
        # saw, and say so in read_basis rather than pretend the number is canonical.
        # These fallbacks already *contain* the chimeric reads, so they are their own
        # denominator — adding the chimeras back would count them twice.
        total_reads = sum(cluster_reads.values()) if cluster_reads else total_size
        read_basis = "cluster" if cluster_reads else "centroid"
        denominator = total_reads
    else:
        # otu_table excludes chimeric reads by construction — chimeric variants are
        # dropped before the mapping step — so they must be added back to ask "what
        # share of the reads that reached clustering was chimeric".
        total_reads = mapped_reads
        read_basis = "otu_table"
        denominator = mapped_reads + chimeric_reads

    return {
        "chimeric_variants": n_chimeric,
        "total_variants": n_total,
        "pct_chimeric_variants": round(100.0 * n_chimeric / n_total, 3) if n_total else 0.0,
        "chimeric_reads": chimeric_reads,
        "total_reads": total_reads,
        "pct_chimeric_reads": round(100.0 * chimeric_reads / denominator, 3) if denominator else 0.0,
        "read_basis": read_basis,
    }


# ── Abundance matrix ─────────────────────────────────────────────────────


def load_matrix(path) -> tuple[list[str], list[str], np.ndarray]:
    """Read a feature x sample abundance table into ``(features, samples, counts)``.

    ``counts`` is shaped ``(n_samples, n_features)`` — the orientation scikit-bio's
    diversity functions expect.
    """
    path = Path(path)
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    with path.open("r", errors="replace") as handle:
        rows = list(csv.reader(handle, delimiter=delimiter))
    if not rows:
        return [], [], np.zeros((0, 0))
    header, *body = rows
    samples = [name.strip() for name in header[1:]]
    features: list[str] = []
    values: list[list[float]] = []
    for row in body:
        if not row or not row[0].strip():
            continue
        features.append(row[0].strip())
        values.append([float(cell or 0) for cell in row[1 : len(samples) + 1]])
    if not features:
        return [], samples, np.zeros((len(samples), 0))
    counts = np.asarray(values, dtype=float).T  # samples x features
    return features, samples, counts


# ── Alpha diversity ──────────────────────────────────────────────────────


def _hill(counts: np.ndarray, q: int) -> float:
    """Hill number of order q: the effective number of species."""
    total = counts.sum()
    if total <= 0:
        return 0.0
    proportions = counts[counts > 0] / total
    if q == 0:
        return float(len(proportions))
    if q == 1:
        return float(math.exp(-np.sum(proportions * np.log(proportions))))
    return float(np.sum(proportions**q) ** (1.0 / (1.0 - q)))


def alpha_diversity_row(sample: str, counts: np.ndarray) -> dict:
    """Alpha diversity of one already-rarefied sample.

    Richness, Chao1, Shannon, Simpson and Pielou come from scikit-bio rather than
    being re-derived here: its Chao1 is the bias-corrected estimator, and its
    Shannon defaults to natural log — the same base the Hill numbers assume, so
    ``hill_q1 == exp(shannon)`` holds exactly.
    """
    from skbio.diversity.alpha import chao1, pielou_e, shannon, simpson

    counts = counts.astype(int)
    total = int(counts.sum())
    observed = int((counts > 0).sum())
    if total == 0:
        return {"sample": sample, "reads": 0, **{column: 0.0 for column in ALPHA_COLUMNS[2:]}}
    return {
        "sample": sample,
        "reads": total,
        "observed": observed,
        "chao1": round(float(chao1(counts)), 3),
        "shannon": round(float(shannon(counts)), 4),
        "simpson": round(float(simpson(counts)), 4),
        "pielou": round(float(pielou_e(counts)) if observed > 1 else 0.0, 4),
        "hill_q0": round(_hill(counts, 0), 3),
        "hill_q1": round(_hill(counts, 1), 3),
        "hill_q2": round(_hill(counts, 2), 3),
    }


def rarefy(counts: np.ndarray, depth: int, rng: np.random.Generator) -> np.ndarray:
    """Subsample one sample's counts down to ``depth`` reads without replacement."""
    total = int(counts.sum())
    if total <= depth:
        return counts.astype(int)
    # Hypergeometric draw over the multiset of reads, feature by feature.
    remaining = total
    left = depth
    out = np.zeros_like(counts, dtype=int)
    for index, value in enumerate(counts.astype(int)):
        if left <= 0:
            break
        drawn = rng.hypergeometric(value, remaining - value, left) if value else 0
        out[index] = drawn
        remaining -= value
        left -= drawn
    return out


def rarefaction_curve(counts: np.ndarray, samples: list[str], rng: np.random.Generator) -> list[dict]:
    """Observed richness against sequencing depth, averaged over random subsamples."""
    rows: list[dict] = []
    for index, sample in enumerate(samples):
        sample_counts = counts[index]
        total = int(sample_counts.sum())
        if total <= 0:
            continue
        depths = np.unique(np.linspace(1, total, num=min(RAREFACTION_POINTS, total)).astype(int))
        for depth in depths:
            richness = [int((rarefy(sample_counts, int(depth), rng) > 0).sum()) for _ in range(RAREFACTION_ITERATIONS)]
            rows.append({"sample": sample, "depth": int(depth), "richness": round(float(np.mean(richness)), 3)})
    return rows


def rank_abundance(counts: np.ndarray, samples: list[str]) -> list[dict]:
    """Whittaker plot data: relative abundance against rank, per sample."""
    rows: list[dict] = []
    for index, sample in enumerate(samples):
        sample_counts = counts[index]
        total = sample_counts.sum()
        if total <= 0:
            continue
        ordered = np.sort(sample_counts[sample_counts > 0])[::-1]
        for rank, value in enumerate(ordered, start=1):
            rows.append({"sample": sample, "rank": rank, "relative_abundance": round(float(value / total), 8)})
    return rows


def feature_intersections(counts: np.ndarray, samples: list[str], max_rows: int = 40) -> list[dict]:
    """UpSet data: how many features occur in exactly each combination of samples.

    An UpSet replaces the Venn diagram because Venn stops being readable past three
    sets, and this scales to dozens of samples. Combinations are stored as a
    ``|``-joined sample list, biggest first.
    """
    if not samples:
        return []
    present = counts > 0
    tally: Counter = Counter()
    for column in range(present.shape[1]):
        members = tuple(sample for index, sample in enumerate(samples) if present[index, column])
        if members:
            tally[members] += 1
    rows = [
        {"combination": "|".join(members), "degree": len(members), "n_features": count}
        for members, count in tally.most_common(max_rows)
    ]
    return rows


OCCUPANCY_COLUMNS = [
    "feature",
    "n_samples_present",
    "samples_present",
    "total_reads",
    "max_sample",
    "max_sample_reads",
    "is_shared",
    "is_unique",
]


def feature_occupancy(counts: np.ndarray, features: list[str], samples: list[str]) -> list[dict]:
    """Which samples each feature was actually found in, and how strongly.

    ``feature_prevalence`` answers "in how many samples?"; this answers "in WHICH
    samples?", which is the question you ask when deciding whether a parasite is
    site-wide or confined to one host. Shared (>1 sample) features sort first, then
    by read total, so the widespread ones are on top.
    """
    if not samples or not features:
        return []
    present = counts > 0
    rows = []
    for index, feature in enumerate(features):
        column = counts[:, index]
        hit_indices = [i for i in range(len(samples)) if present[i, index]]
        if not hit_indices:
            continue
        top = max(hit_indices, key=lambda i: column[i])
        rows.append(
            {
                "feature": feature,
                "n_samples_present": len(hit_indices),
                "samples_present": ";".join(samples[i] for i in hit_indices),
                "total_reads": int(column.sum()),
                "max_sample": samples[top],
                "max_sample_reads": int(column[top]),
                "is_shared": int(len(hit_indices) > 1),
                "is_unique": int(len(hit_indices) == 1),
            }
        )
    rows.sort(key=lambda row: (-row["n_samples_present"], -row["total_reads"]))
    return rows


def feature_prevalence(counts: np.ndarray, features: list[str], core_fraction: float = DEFAULT_CORE_FRACTION) -> list[dict]:
    """Per-feature prevalence and whether it clears the core threshold."""
    n_samples = counts.shape[0]
    if n_samples == 0:
        return []
    present = counts > 0
    rows = []
    for index, feature in enumerate(features):
        n_present = int(present[:, index].sum())
        rows.append(
            {
                "feature": feature,
                "n_samples_present": n_present,
                "prevalence": round(n_present / n_samples, 4),
                "total_reads": int(counts[:, index].sum()),
                "is_core": int(n_present / n_samples >= core_fraction),
            }
        )
    return rows


# ── Persistence ──────────────────────────────────────────────────────────


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def build_analysis_metrics(
    analysis_dir,
    matrix_path,
    *,
    otus_fasta: Optional[Path] = None,
    chimeras_fasta: Optional[Path] = None,
    otus_raw_fasta: Optional[Path] = None,
    clusters_uc: Optional[Path] = None,
    core_fraction: float = DEFAULT_CORE_FRACTION,
    seed: int = 0,
    metrics_dir: Optional[Path] = None,
) -> dict[str, Path]:
    """Compute and persist every analysis-stage table. Returns name -> path.

    ``metrics_dir`` lets a caller keep each feature set's tables in its own
    directory (``metrics/otu``, ``metrics/asv``). OTUs and ASVs are different
    units and must never land in the same table.
    """
    analysis_dir = Path(analysis_dir)
    metrics_dir = Path(metrics_dir) if metrics_dir else analysis_dir / "metrics"
    written: dict[str, Path] = {}

    if chimeras_fasta and otus_raw_fasta:
        # matrix_path is what makes total_reads the canonical, mapping-based total
        # rather than the clustering's own smaller population.
        stats = chimera_stats(chimeras_fasta, otus_raw_fasta, clusters_uc, matrix_path)
        written["chimera_stats_csv"] = _write_csv(metrics_dir / "chimera_stats.csv", list(stats), [stats])

    if otus_fasta and Path(otus_fasta).exists():
        rows = [{"length": length, "count": count} for length, count in sorted(fasta_lengths(otus_fasta).items())]
        if rows:
            written["feature_length_csv"] = _write_csv(metrics_dir / "feature_length_distribution.csv", ["length", "count"], rows)

    if not (matrix_path and Path(matrix_path).exists()):
        return written

    features, samples, counts = load_matrix(matrix_path)
    if not features or not samples:
        return written

    rng = np.random.default_rng(seed)
    depths = counts.sum(axis=1)
    non_empty = depths > 0
    # Rarefy to the shallowest non-empty sample: every comparison is then made at
    # equal sequencing effort, which is the only way alpha/beta are comparable.
    depth = int(depths[non_empty].min()) if non_empty.any() else 0

    written["rarefaction_depth_csv"] = _write_csv(
        metrics_dir / "rarefaction_depth.csv",
        ["depth", "n_samples", "n_samples_dropped"],
        [{"depth": depth, "n_samples": int(non_empty.sum()), "n_samples_dropped": int((~non_empty).sum())}],
    )

    alpha_rows = []
    for index, sample in enumerate(samples):
        if not non_empty[index]:
            continue
        alpha_rows.append(alpha_diversity_row(sample, rarefy(counts[index], depth, rng)))
    if alpha_rows:
        written["alpha_diversity_csv"] = _write_csv(metrics_dir / "alpha_diversity.csv", ALPHA_COLUMNS, alpha_rows)

    curve = rarefaction_curve(counts, samples, rng)
    if curve:
        written["rarefaction_csv"] = _write_csv(metrics_dir / "rarefaction.csv", ["sample", "depth", "richness"], curve)

    ranks = rank_abundance(counts, samples)
    if ranks:
        written["rank_abundance_csv"] = _write_csv(
            metrics_dir / "rank_abundance.csv", ["sample", "rank", "relative_abundance"], ranks
        )

    prevalence = feature_prevalence(counts, features, core_fraction)
    if prevalence:
        written["prevalence_csv"] = _write_csv(
            metrics_dir / "feature_prevalence.csv",
            ["feature", "n_samples_present", "prevalence", "total_reads", "is_core"],
            prevalence,
        )

    intersections = feature_intersections(counts, samples)
    if intersections:
        written["intersections_csv"] = _write_csv(
            metrics_dir / "feature_intersections.csv", ["combination", "degree", "n_features"], intersections
        )

    occupancy = feature_occupancy(counts, features, samples)
    if occupancy:
        written["occupancy_csv"] = _write_csv(metrics_dir / "feature_occupancy.csv", OCCUPANCY_COLUMNS, occupancy)

    return written


def _path_if_exists(value) -> Optional[Path]:
    return Path(value) if value and Path(str(value)).exists() else None


def feature_sets_from_outputs(analysis_outputs: dict) -> list[dict]:
    """Describe every feature set the analysis stage actually produced.

    An OTU and an ASV are different units of biology — a 97 %-identity cluster
    versus an exact denoised sequence — so each gets its own entry here, and every
    downstream table is built per entry. They are never concatenated.
    """
    sets: list[dict] = []

    otu_fasta = _path_if_exists(analysis_outputs.get("combined_fasta_vsearch"))
    otu_matrix = _path_if_exists(analysis_outputs.get("otu_matrix") or analysis_outputs.get("otu_table_tsv"))
    if otu_fasta:
        sets.append({"key": "otu", "label": FEATURE_SET_LABELS["otu"], "fasta": otu_fasta, "matrix": otu_matrix})

    dada2_files = (analysis_outputs.get("dada2") or {}).get("files") or {}
    asv_fasta = _path_if_exists(analysis_outputs.get("combined_fasta_dada2") or dada2_files.get("asv_sequences_fasta"))
    asv_matrix = _path_if_exists(dada2_files.get("asv_abundance_csv"))
    if asv_fasta:
        sets.append({"key": "asv", "label": FEATURE_SET_LABELS["asv"], "fasta": asv_fasta, "matrix": asv_matrix})

    return sets


def build_metrics_from_outputs(analysis_dir, analysis_outputs: dict) -> dict[str, str]:
    """Persist the analysis tables for every feature set the run produced.

    Returns a flat name -> path map. Per-feature-set keys are namespaced
    (``otu.alpha_diversity_csv``, ``asv.alpha_diversity_csv``); the unprefixed keys
    are kept as aliases for the primary feature set so existing readers still work.
    """
    pool_files = (analysis_outputs.get("vsearch_pool") or {}).get("files") or {}
    analysis_dir = Path(analysis_dir)
    written: dict[str, str] = {}

    for index, feature_set in enumerate(feature_sets_from_outputs(analysis_outputs)):
        key = feature_set["key"]
        # Chimera accounting only exists for the VSEARCH branch; DADA2 removes
        # bimeras inside its own run and reports them in read_tracking.csv.
        is_otu = key == "otu"
        tables = build_analysis_metrics(
            analysis_dir,
            feature_set["matrix"],
            otus_fasta=feature_set["fasta"],
            chimeras_fasta=_path_if_exists(pool_files.get("chimeras_fasta")) if is_otu else None,
            otus_raw_fasta=_path_if_exists(pool_files.get("otus_raw_fasta")) if is_otu else None,
            clusters_uc=_path_if_exists(pool_files.get("clusters_uc")) if is_otu else None,
            metrics_dir=analysis_dir / "metrics" / key,
        )
        for name, path in tables.items():
            written[f"{key}.{name}"] = str(path)
            # The first feature set is the primary one; alias it so callers that
            # predate feature sets keep resolving to a real table.
            if index == 0:
                written.setdefault(name, str(path))

    return written


def analysis_metrics_paths(analysis_dir, feature_set: str = "otu") -> dict[str, Path]:
    """Return the persisted analysis tables for one feature set."""
    analysis_dir = Path(analysis_dir)
    names = {
        "chimera_stats_csv": "chimera_stats.csv",
        "feature_length_csv": "feature_length_distribution.csv",
        "rarefaction_depth_csv": "rarefaction_depth.csv",
        "alpha_diversity_csv": "alpha_diversity.csv",
        "rarefaction_csv": "rarefaction.csv",
        "rank_abundance_csv": "rank_abundance.csv",
        "prevalence_csv": "feature_prevalence.csv",
        "intersections_csv": "feature_intersections.csv",
        "occupancy_csv": "feature_occupancy.csv",
    }
    # Runs predating feature sets wrote straight into metrics/.
    for metrics_dir in (analysis_dir / "metrics" / feature_set, analysis_dir / "metrics"):
        found = {key: metrics_dir / name for key, name in names.items() if (metrics_dir / name).exists()}
        if found:
            return found
    return {}
