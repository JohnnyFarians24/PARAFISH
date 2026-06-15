"""Parsers for FastQC output artifacts."""

from __future__ import annotations

import os
import zipfile
from pathlib import Path


def _parse_fastqc_basic_stats_lines(lines):
    stats = {}
    in_basic = False
    for raw_line in lines:
        line = raw_line.strip()
        if line.startswith(">>Basic Statistics"):
            in_basic = True
            continue
        if in_basic and line.startswith("#"):
            continue
        if in_basic and line.startswith(">>END_MODULE"):
            break
        if in_basic and "\t" in line:
            key, value = line.split("\t", 1)
            stats[key.strip()] = value.strip()
    return stats


def _parse_fastqc_basic_stats_from_zip(zip_path):
    if not zip_path or not os.path.exists(zip_path) or os.path.getsize(zip_path) == 0:
        return {}
    try:
        with zipfile.ZipFile(zip_path, "r") as archive:
            for member in archive.namelist():
                if member.endswith("fastqc_data.txt"):
                    with archive.open(member) as handle:
                        lines = [line.decode("utf-8", errors="replace") for line in handle.readlines()]
                    return _parse_fastqc_basic_stats_lines(lines)
    except Exception:
        return {}
    return {}


def parse_fastqc_report(report_path: Path) -> dict[str, str]:
    """Parse a FastQC zip report into a normalized summary dictionary."""
    report_path = Path(report_path)

    if not report_path.exists():
        return {"report_path": str(report_path), "error": "File not found"}

    if report_path.suffix == ".html":
        zip_path = Path(str(report_path).replace("_fastqc.html", "_fastqc.zip"))
    else:
        zip_path = report_path

    stats = _parse_fastqc_basic_stats_from_zip(zip_path)
    return {
        "report_path": str(report_path),
        "zip_path": str(zip_path),
        **stats,
    }


def collect_fastqc_reports(output_dir: str | Path) -> dict[str, dict]:
    """
    Walk a QC output directory and return parsed stats per sample.

    Returns:
        {"sample_id": {"r1": {...stats}, "r2": {...stats}}}
    """
    output_dir = Path(output_dir)
    if not output_dir.exists():
        return {}

    results: dict[str, dict] = {}

    for zip_file in sorted(output_dir.rglob("*_fastqc.zip")):
        name = zip_file.stem
        stats = parse_fastqc_report(zip_file)

        if "_R1_fastqc" in name:
            sample_id = name.replace("_R1_fastqc", "")
            strand = "r1"
        elif "_R2_fastqc" in name:
            sample_id = name.replace("_R2_fastqc", "")
            strand = "r2"
        else:
            sample_id = name.replace("_fastqc", "")
            strand = "r1"

        results.setdefault(sample_id, {})
        results[sample_id][strand] = stats

    return results
