"""Parsers for FastQC output artifacts."""

from __future__ import annotations

import os
import zipfile
from typing import Optional


def _read_fastqc_data(zip_path) -> str:
    """Return the raw text of ``fastqc_data.txt`` inside a FastQC zip, or ""."""
    if not zip_path or not os.path.exists(zip_path) or os.path.getsize(zip_path) == 0:
        return ""
    try:
        with zipfile.ZipFile(zip_path, "r") as archive:
            members = [name for name in archive.namelist() if name.endswith("fastqc_data.txt")]
            if not members:
                return ""
            return archive.read(members[0]).decode("utf-8", errors="replace")
    except (zipfile.BadZipFile, KeyError, OSError):
        return ""


def _split_modules(text: str) -> dict[str, tuple[str, list[list[str]]]]:
    """Split ``fastqc_data.txt`` into ``{module: (status, rows)}``.

    Each module is delimited by ``>>Name<TAB>status`` ... ``>>END_MODULE``. Rows are
    the tab-split data lines; ``#`` header lines are dropped, except that FastQC
    smuggles ``#Total Deduplicated Percentage`` in as a comment, so comment lines
    that carry a value are kept under a leading ``#`` marker.
    """
    modules: dict[str, tuple[str, list[list[str]]]] = {}
    name: Optional[str] = None
    status = ""
    rows: list[list[str]] = []
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith(">>END_MODULE"):
            if name is not None:
                modules[name] = (status, rows)
            name, status, rows = None, "", []
            continue
        if line.startswith(">>"):
            header = line[2:].split("\t")
            name = header[0].strip()
            status = header[1].strip().lower() if len(header) > 1 else ""
            rows = []
            continue
        if name is None or not line.strip():
            continue
        if line.startswith("#"):
            fields = line.lstrip("#").split("\t")
            if len(fields) > 1:
                rows.append(["#" + fields[0].strip(), *[f.strip() for f in fields[1:]]])
            continue
        rows.append([field.strip() for field in line.split("\t")])
    if name is not None:
        modules[name] = (status, rows)
    return modules


def _to_float(value) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _midpoint(label: str) -> Optional[float]:
    """FastQC bins positions and lengths as "10" or "150-151"; take the midpoint."""
    text = str(label).strip()
    if "-" in text:
        low, _, high = text.partition("-")
        low_value, high_value = _to_float(low), _to_float(high)
        if low_value is None or high_value is None:
            return None
        return (low_value + high_value) / 2.0
    return _to_float(text)


def parse_fastqc_data(zip_path) -> dict:
    """Parse the full FastQC data file into a structured payload.

    Returns ``{}`` when the zip is missing or unreadable. Keys:

    ``module_status``      {module name: pass|warn|fail}
    ``total_sequences``    read count
    ``pct_gc``             overall %GC
    ``mean_length``        abundance-weighted mean read length
    ``pct_duplication``    100 - "Total Deduplicated Percentage"
    ``pct_adapter_max``    highest adapter content reached at any position
    ``pct_reads_q30``      share of reads whose *mean* quality is >= 30
    ``per_base_quality``   [(position, mean_quality)]
    ``length_distribution``[(length, count)]
    """
    text = _read_fastqc_data(zip_path)
    if not text:
        return {}
    modules = _split_modules(text)
    payload: dict = {"module_status": {name: status for name, (status, _) in modules.items()}}

    basic_rows = modules.get("Basic Statistics", ("", []))[1]
    basic = {row[0]: row[1] for row in basic_rows if len(row) > 1}
    payload["total_sequences"] = _to_float(basic.get("Total Sequences"))
    payload["pct_gc"] = _to_float(basic.get("%GC"))

    per_base = []
    for row in modules.get("Per base sequence quality", ("", []))[1]:
        if len(row) > 1 and not row[0].startswith("#"):
            position, mean = _midpoint(row[0]), _to_float(row[1])
            if position is not None and mean is not None:
                per_base.append((position, mean))
    payload["per_base_quality"] = per_base

    # FastQC reports per-read *mean* quality, so this is the share of reads whose
    # average quality clears Q30 — not the share of bases. Naming it anything else
    # would misrepresent what FastQC measured.
    q_counts = []
    for row in modules.get("Per sequence quality scores", ("", []))[1]:
        if len(row) > 1 and not row[0].startswith("#"):
            quality, count = _to_float(row[0]), _to_float(row[1])
            if quality is not None and count is not None:
                q_counts.append((quality, count))
    total_q = sum(count for _, count in q_counts)
    payload["pct_reads_q30"] = (
        round(100.0 * sum(count for quality, count in q_counts if quality >= 30) / total_q, 2) if total_q else None
    )

    lengths = []
    for row in modules.get("Sequence Length Distribution", ("", []))[1]:
        if len(row) > 1 and not row[0].startswith("#"):
            length, count = _midpoint(row[0]), _to_float(row[1])
            if length is not None and count is not None:
                lengths.append((length, count))
    payload["length_distribution"] = lengths
    total_reads = sum(count for _, count in lengths)
    payload["mean_length"] = (
        round(sum(length * count for length, count in lengths) / total_reads, 1) if total_reads else None
    )

    dedup = None
    for row in modules.get("Sequence Duplication Levels", ("", []))[1]:
        if row[0] == "#Total Deduplicated Percentage" and len(row) > 1:
            dedup = _to_float(row[1])
            break
    payload["pct_duplication"] = round(100.0 - dedup, 2) if dedup is not None else None

    adapter_max = 0.0
    for row in modules.get("Adapter Content", ("", []))[1]:
        if row[0].startswith("#"):
            continue
        for cell in row[1:]:
            value = _to_float(cell)
            if value is not None:
                adapter_max = max(adapter_max, value)
    payload["pct_adapter_max"] = round(adapter_max, 2) if modules.get("Adapter Content") else None

    return payload
