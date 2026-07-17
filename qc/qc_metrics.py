"""Persisted QC metrics.

The QC stage writes its tables once, here, straight after FastQC runs. The Dash
callbacks only read these CSVs back — they never re-parse the FastQC zips, so
switching pages does not re-do the work.
"""

from __future__ import annotations

import csv
from pathlib import Path
from statistics import median
from typing import Optional

from qc.fastqc_parser import parse_fastqc_data

# A sample whose depth falls below this fraction of the run's median depth is
# flagged: it will dominate the rarefaction floor and drag every alpha/beta
# comparison down with it, so it is a candidate for removal before diversity.
DEPTH_OUTLIER_FRACTION = 0.5

OVERVIEW_COLUMNS = [
    "sample",
    "read",
    "total_reads",
    "pct_gc",
    "mean_length",
    "pct_duplication",
    "pct_adapter_max",
    "pct_reads_q30",
    "n_fail",
    "n_warn",
    "is_depth_outlier",
]


def _strand_of(zip_path: Path) -> str:
    name = zip_path.stem
    if "_R2_fastqc" in name:
        return "R2"
    return "R1"


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _collect(qc_dir: Path) -> list[tuple[str, str, dict]]:
    """Return ``(sample_id, strand, payload)`` for every FastQC zip under ``qc_dir``.

    The sample id comes from the containing directory, which the runner names after
    ``Sample.internal_id``. Deriving it from the FASTQ filename instead would report
    "APTL64B" where the rest of the pipeline says "aptl64b_02".
    """
    collected: list[tuple[str, str, dict]] = []
    for zip_path in sorted(qc_dir.rglob("*_fastqc.zip")):
        payload = parse_fastqc_data(zip_path)
        if not payload:
            continue
        collected.append((zip_path.parent.name, _strand_of(zip_path), payload))
    return collected


def _depth_outliers(entries: list[tuple[str, str, dict]]) -> set[str]:
    depths = {
        sample: payload.get("total_sequences")
        for sample, strand, payload in entries
        if strand == "R1" and payload.get("total_sequences")
    }
    if len(depths) < 2:
        return set()
    threshold = DEPTH_OUTLIER_FRACTION * median(depths.values())
    return {sample for sample, depth in depths.items() if depth < threshold}


def build_qc_metrics(qc_dir) -> dict[str, Path]:
    """Parse every FastQC report under ``qc_dir`` and persist the QC tables.

    Returns a mapping of logical name -> written path (empty when no reports).
    """
    qc_dir = Path(qc_dir)
    if not qc_dir.exists():
        return {}
    entries = _collect(qc_dir)
    if not entries:
        return {}

    outliers = _depth_outliers(entries)
    metrics_dir = qc_dir / "metrics"

    overview: list[dict] = []
    modules: list[dict] = []
    per_base: list[dict] = []
    lengths: list[dict] = []

    for sample, strand, payload in entries:
        status = payload.get("module_status") or {}
        overview.append(
            {
                "sample": sample,
                "read": strand,
                "total_reads": int(payload["total_sequences"]) if payload.get("total_sequences") else "",
                "pct_gc": payload.get("pct_gc") if payload.get("pct_gc") is not None else "",
                "mean_length": payload.get("mean_length") if payload.get("mean_length") is not None else "",
                "pct_duplication": payload.get("pct_duplication") if payload.get("pct_duplication") is not None else "",
                "pct_adapter_max": payload.get("pct_adapter_max") if payload.get("pct_adapter_max") is not None else "",
                "pct_reads_q30": payload.get("pct_reads_q30") if payload.get("pct_reads_q30") is not None else "",
                "n_fail": sum(1 for value in status.values() if value == "fail"),
                "n_warn": sum(1 for value in status.values() if value == "warn"),
                "is_depth_outlier": int(sample in outliers),
            }
        )
        for module, value in status.items():
            modules.append({"sample": sample, "read": strand, "module": module, "status": value})
        for position, quality in payload.get("per_base_quality") or []:
            per_base.append({"sample": sample, "read": strand, "position": position, "mean_quality": round(quality, 3)})
        for length, count in payload.get("length_distribution") or []:
            lengths.append({"sample": sample, "read": strand, "length": length, "count": int(count)})

    return {
        "qc_overview_csv": _write_csv(metrics_dir / "qc_overview.csv", OVERVIEW_COLUMNS, overview),
        "qc_module_status_csv": _write_csv(metrics_dir / "qc_module_status.csv", ["sample", "read", "module", "status"], modules),
        "qc_per_base_quality_csv": _write_csv(
            metrics_dir / "qc_per_base_quality.csv", ["sample", "read", "position", "mean_quality"], per_base
        ),
        "qc_length_distribution_csv": _write_csv(
            metrics_dir / "qc_length_distribution.csv", ["sample", "read", "length", "count"], lengths
        ),
    }


def qc_metrics_paths(qc_dir) -> dict[str, Path]:
    """Return the persisted QC tables that exist under ``qc_dir``."""
    metrics_dir = Path(qc_dir) / "metrics"
    names = {
        "qc_overview_csv": "qc_overview.csv",
        "qc_module_status_csv": "qc_module_status.csv",
        "qc_per_base_quality_csv": "qc_per_base_quality.csv",
        "qc_length_distribution_csv": "qc_length_distribution.csv",
    }
    return {key: metrics_dir / name for key, name in names.items() if (metrics_dir / name).exists()}
