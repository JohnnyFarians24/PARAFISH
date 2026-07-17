"""Cross-stage results integration.

The Results page answers questions no single stage can: where did the reads go
from raw FASTQ to classified feature, how do samples compare on one line each,
and — for a parasitologist — where does each target taxon turn up and how much
to trust it.

Everything here reads the per-stage ``metrics/`` CSVs that stages 1-5 already
wrote, so nothing is recomputed; the functions just join them.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional


def _read_rows(path) -> list[dict]:
    path = Path(path) if path else None
    if not path or not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle))


def _to_int(value) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return 0


# Every stage writes one directory per feature set (``databases/otu/``,
# ``analysis/metrics/otu/``). The resolvers in ``analysis_metrics`` / ``database_metrics``
# own that layout, including the fallback to the flat pre-feature-set one. This module
# must never re-derive it: re-deriving it, and missing the feature-set level, is exactly
# what made every number on the Results page read blank on a run that had classified
# 9/9 features (DEVLOG 2026-07-17).
#
# Like the rest of the Results layer, this reads the OTU set; a DADA2-only run is a
# known pre-existing gap, not one introduced here.
DEFAULT_FEATURE_SET = "otu"


def _analysis_metric(run_output_dir: Path, key: str) -> Optional[Path]:
    from analysis.analysis_metrics import analysis_metrics_paths

    return analysis_metrics_paths(run_output_dir / "analysis", DEFAULT_FEATURE_SET).get(key)


def _database_metric(run_output_dir: Path, key: str) -> Optional[Path]:
    from databases.database_metrics import database_metrics_paths

    return database_metrics_paths(run_output_dir / "databases", DEFAULT_FEATURE_SET).get(key)


def _taxonomy_table(run_output_dir: Path, filename: str) -> Optional[Path]:
    from databases.database_metrics import taxonomy_table_path

    return taxonomy_table_path(run_output_dir / "databases", filename, DEFAULT_FEATURE_SET)


# ── End-to-end read fate (Sankey) ────────────────────────────────────────

# Ordered pipeline stops a read passes through. Each stop's value is summed across
# samples; the flows between consecutive stops are what the Sankey draws.
_FATE_STAGES = ["Raw", "Trimmed", "Merged", "Filtered", "Clustered", "Classified"]


def read_fate_sankey(run_output_dir, *, feature_reads: Optional[int] = None, classified_reads: Optional[int] = None) -> dict:
    """Whole-pipeline read counts, for a Sankey from raw reads to classified reads.

    Reads come from the preprocess ``read_tracking.csv`` (raw → trimmed → merged);
    the later stops are supplied by the caller from the analysis/database metrics,
    since only they know the clustered and classified totals.
    """
    run_output_dir = Path(run_output_dir) if run_output_dir else None
    if run_output_dir is None:
        return {}
    tracking = _read_rows(run_output_dir / "preprocess" / "metrics" / "read_tracking.csv")
    if not tracking:
        return {}

    raw = sum(_to_int(row.get("input_reads")) for row in tracking)
    trimmed = sum(_to_int(row.get("trimmed_reads")) for row in tracking)
    merged = sum(_to_int(row.get("assembled_reads") or row.get("surviving_reads")) for row in tracking)

    stages = ["Raw", "Trimmed", "Merged"]
    values = [raw, trimmed, merged]
    if feature_reads is not None:
        stages.append("Clustered")
        values.append(int(feature_reads))
    if classified_reads is not None:
        stages.append("Classified")
        values.append(int(classified_reads))

    # Consecutive flows, never letting a downstream stop exceed its upstream (a
    # merged count above the trimmed count would mean the accounting is wrong).
    flows = []
    for index in range(len(stages) - 1):
        source_value = values[index]
        target_value = min(values[index + 1], source_value) if source_value else values[index + 1]
        flows.append({"source": stages[index], "target": stages[index + 1], "value": target_value})
    return {"stages": stages, "values": values, "flows": flows}


def clustered_read_total(run_output_dir) -> Optional[int]:
    """Reads that survived into the clustered feature set.

    Read from the chimera table the analysis stage persisted. Recomputing it here
    would mean re-scanning ``clusters.uc``, which holds one row per dereplicated
    read — fine on a toy run, minutes on a real one.
    """
    run_output_dir = Path(run_output_dir) if run_output_dir else None
    if run_output_dir is None:
        return None
    rows = _read_rows(_analysis_metric(run_output_dir, "chimera_stats_csv"))
    if not rows:
        return None
    total = _to_int(rows[0].get("total_reads"))
    chimeric = _to_int(rows[0].get("chimeric_reads"))
    return max(total - chimeric, 0) or None


def species_classified_total(run_output_dir) -> Optional[int]:
    """Reads on features that reached a species-level consensus.

    Read from the taxonomic-resolution table the databases stage persisted.
    """
    run_output_dir = Path(run_output_dir) if run_output_dir else None
    if run_output_dir is None:
        return None
    for row in _read_rows(_database_metric(run_output_dir, "resolution_csv")):
        if row.get("rank") == "species":
            return _to_int(row.get("n_reads"))
    return None


# ── Per-sample dashboard ─────────────────────────────────────────────────


def sample_dashboard(
    run_output_dir,
    control_samples: Optional[list[str]] = None,
    target_taxa: Optional[list[str]] = None,
) -> list[dict]:
    """One row per sample: reads, richness, diversity, % classified, target hits.

    Joins the preprocess waterfall, the analysis alpha-diversity table and the
    taxonomy consensus. Control samples are marked, never silently dropped.
    """
    run_output_dir = Path(run_output_dir) if run_output_dir else None
    if run_output_dir is None:
        return []
    controls = set(control_samples or [])

    surviving = {row["sample"]: _to_int(row.get("surviving_reads")) for row in _read_rows(run_output_dir / "preprocess" / "metrics" / "read_tracking.csv")}
    alpha = {row["sample"]: row for row in _read_rows(_analysis_metric(run_output_dir, "alpha_diversity_csv"))}

    pct_classified, target_hits = _classification_by_sample(run_output_dir, target_taxa)

    samples = sorted(set(surviving) | set(alpha) | set(pct_classified))
    rows = []
    for sample in samples:
        entry = alpha.get(sample, {})
        rows.append(
            {
                "sample": sample,
                "is_control": int(sample in controls),
                "reads": surviving.get(sample, _to_int(entry.get("reads"))),
                "richness": _to_int(entry.get("observed")),
                "shannon": float(entry["shannon"]) if entry.get("shannon") else 0.0,
                "pct_species": pct_classified.get(sample, 0.0),
                "n_targets": len(target_hits.get(sample, set())),
            }
        )
    return rows


def _classification_by_sample(run_output_dir: Path, target_taxa: Optional[list[str]]) -> tuple[dict, dict]:
    """Per-sample % of reads classified to species and which target taxa appear.

    Reads are attributed to a sample through the OTU abundance matrix, then to a
    taxon through the consensus. Without both, this degrades to empty rather than
    guessing.
    """
    from analysis.analysis_metrics import load_matrix
    from databases.database_metrics import _CONSENSUS_COLUMN, _blank, _read_consensus

    matrix_path = run_output_dir / "analysis" / "otu_matrix.csv"
    consensus_path = _taxonomy_table(run_output_dir, "taxonomy_consensus.csv")
    if not matrix_path.exists() or consensus_path is None:
        return {}, {}

    consensus = _read_consensus(consensus_path)
    features, samples, counts = load_matrix(matrix_path)
    if not features:
        return {}, {}

    targets = {taxon.strip().lower() for taxon in (target_taxa or []) if taxon.strip()}
    species_column = _CONSENSUS_COLUMN["species"]

    pct_classified: dict[str, float] = {}
    target_hits: dict[str, set] = {}
    # counts is (samples x features); axis=1 gives the reads *per sample*. Summing
    # over axis=0 would give per-feature totals and make the percentage meaningless
    # (it read as 650% before this was fixed).
    totals = counts.sum(axis=1)
    for sample_index, sample in enumerate(samples):
        total = totals[sample_index]
        classified = 0.0
        for feature_index, feature in enumerate(features):
            reads = counts[sample_index, feature_index]
            if reads <= 0:
                continue
            row = consensus.get(feature, {})
            if not _blank(row.get(species_column)):
                classified += reads
            if targets:
                for column in _CONSENSUS_COLUMN.values():
                    label = str(row.get(column, "") or "").strip().lower()
                    if label in targets:
                        target_hits.setdefault(sample, set()).add(label)
        pct_classified[sample] = round(100.0 * float(classified) / float(total), 2) if total else 0.0
    return pct_classified, target_hits


# ── Target-taxon report card ─────────────────────────────────────────────


def target_taxon_report(run_output_dir, target_taxa: list[str], control_samples: Optional[list[str]] = None) -> list[dict]:
    """For each requested taxon: where it appears, how prevalent, how well supported.

    Prevalence (with its Wilson interval) and reliability are read straight from
    the taxonomy stage's tables so the numbers match what the Databases page shows.
    """
    run_output_dir = Path(run_output_dir) if run_output_dir else None
    targets = [taxon.strip() for taxon in (target_taxa or []) if taxon.strip()]
    if run_output_dir is None or not targets:
        return []

    prevalence = {row["taxon"].strip().lower(): row for row in _read_rows(_database_metric(run_output_dir, "taxon_prevalence_csv"))}
    _, per_sample_hits = _classification_by_sample(run_output_dir, targets)
    reliability = _reliability_by_taxon(run_output_dir)

    rows = []
    for taxon in targets:
        key = taxon.lower()
        prev = prevalence.get(key, {})
        best_reliability = _best_reliability_for(key, reliability)
        samples_with = sorted(sample for sample, hits in per_sample_hits.items() if key in hits and sample not in set(control_samples or []))
        rows.append(
            {
                "taxon": taxon,
                "detected_in": "; ".join(samples_with) if samples_with else "—",
                "n_samples_detected": _to_int(prev.get("n_detected")) or len(samples_with),
                "n_samples": _to_int(prev.get("n_samples")),
                "prevalence": float(prev["prevalence"]) if prev.get("prevalence") else None,
                "ci_low": float(prev["ci_low"]) if prev.get("ci_low") else None,
                "ci_high": float(prev["ci_high"]) if prev.get("ci_high") else None,
                "total_reads": _to_int(prev.get("total_reads")),
                "best_reliability": best_reliability,
            }
        )
    return rows


def _best_reliability_for(target_key: str, reliability: dict[str, str]) -> str:
    """Best grade among assignments whose taxon is, or contains, the target.

    A genus target (``diplostomum``) must pick up its species assignments
    (``diplostomum spathaceum``); an exact-key lookup would miss them.
    """
    best = "—"
    for taxon, grade in reliability.items():
        if taxon == target_key or target_key in taxon.split():
            if best == "—" or _RELIABILITY_ORDER[grade] > _RELIABILITY_ORDER.get(best, -1):
                best = grade
    return best


_RELIABILITY_ORDER = {"high": 3, "moderate": 2, "low": 1, "unknown": 0, "unassigned": 0}


def _reliability_by_taxon(run_output_dir: Path) -> dict[str, str]:
    """Best reliability grade seen for each taxon across all source assignments."""
    best: dict[str, str] = {}
    for row in _read_rows(_taxonomy_table(run_output_dir, "source_assignments.csv")):
        taxon = str(row.get("taxon", "") or "").strip().lower()
        grade = str(row.get("reliability", "") or "").strip().lower()
        if not taxon or grade not in _RELIABILITY_ORDER:
            continue
        if taxon not in best or _RELIABILITY_ORDER[grade] > _RELIABILITY_ORDER[best[taxon]]:
            best[taxon] = grade
    return best
