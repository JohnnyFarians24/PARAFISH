"""Sample naming parsing and manifest row validation for analysis inputs."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Optional, Tuple

from core.config import SUPPORTED_MARKERS, SUPPORTED_MODES


@dataclass
class ParsedFastqName:
    """Normalized sample naming data used across analysis workflows."""

    display_label: str
    internal_sample_id: str
    marker: Optional[str]
    mode: str
    mate: Optional[int] = None


def _infer_marker_from_name(name: str) -> Optional[str]:
    lowered = str(name or "").lower()
    if any(token in lowered for token in ("coi", "cox1", "mtcoi")):
        return "COI"
    if any(token in lowered for token in ("16s", "16s_r", "16s_rna")):
        return "16S"
    if any(token in lowered for token in ("18s", "18s_r", "euk")):
        return "18S"
    if "its" in lowered:
        return "ITS"
    if re.search(r"(?i)(?:^|[^A-Za-z0-9])I\d+[A-Z]?(?:$|[^A-Za-z0-9])", name) and not any(
        token in lowered for token in ("coi", "cox1", "its", "16s", "28s")
    ):
        return "18S"
    return None


def _strip_extensions(filename: str) -> str:
    stem = os.path.basename(str(filename or "")).strip()
    for ext in (
        ".fastq.gz",
        ".fq.gz",
        ".fastq",
        ".fq",
        ".fasta.gz",
        ".fa.gz",
        ".fasta",
        ".fa",
        ".gz",
    ):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    return stem


def _detect_mode_and_mate(stem: str) -> Tuple[str, Optional[int]]:
    lowered = stem.lower()
    if "extendedfrag" in lowered or "assembled" in lowered or "merged" in lowered:
        return "SE", None
    paired_patterns = [
        (r"(?:[_\.-]r?1|[_\.-]forward|[_\.-]fwd)$", 1),
        (r"(?:[_\.-]r?2|[_\.-]reverse|[_\.-]rev)$", 2),
    ]
    for pattern, mate in paired_patterns:
        if re.search(pattern, lowered):
            return "PE", mate
    return "SE", None


def parse_fastq_name(raw_path: str) -> ParsedFastqName:
    """Parse a FASTQ/FASTA file name into display and internal naming fields."""
    stem = _strip_extensions(raw_path)
    mode, mate = _detect_mode_and_mate(stem)
    normalized = re.sub(r"(?i)(?:[_\.-]?(?:r|fwd|rev)?[_\.-]?(?:1|2|forward|reverse))$", "", stem)
    normalized = normalized.strip("_.-") or stem
    marker = _infer_marker_from_name(stem)
    internal_sample_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", normalized).strip("._-") or "sample"
    display_label = normalized.replace("_", " ").strip() or internal_sample_id
    return ParsedFastqName(
        display_label=display_label,
        internal_sample_id=internal_sample_id,
        marker=marker,
        mode=mode,
        mate=mate,
    )


def validate_manifest_rows(manifest_rows: list[dict[str, str]], analysis_mode: str) -> tuple[bool, str]:
    """Validate the analysis manifest rows expected by the current workflow."""
    if not isinstance(manifest_rows, list) or not manifest_rows:
        return False, "Manifest is empty."

    allowed_modes = set(SUPPORTED_MODES)
    required_keys = {"internal_sample_id", "display_label", "mode"}
    for index, row in enumerate(manifest_rows, start=1):
        if not isinstance(row, dict):
            return False, f"Manifest row {index} is invalid."
        missing = required_keys.difference(row.keys())
        if missing:
            return False, f"Manifest row {index} is missing keys: {', '.join(sorted(missing))}."
        mode = str(row.get("mode") or "SE").upper()
        if mode not in allowed_modes:
            return False, f"Manifest row {index} has unsupported mode '{mode}'."
        if not str(row.get("internal_sample_id") or "").strip():
            return False, f"Manifest row {index} is missing an internal_sample_id."
        if not str(row.get("display_label") or "").strip():
            return False, f"Manifest row {index} is missing a display_label."
        if mode == "PE" and not (row.get("forward") and row.get("reverse")):
            return False, f"Manifest row {index} requires both forward and reverse reads for paired-end mode."
        if mode == "SE" and not (row.get("single") or row.get("forward")):
            return False, f"Manifest row {index} requires a single-end input."

    mode_value = str(analysis_mode or "").strip().lower()
    if mode_value not in {"vsearch_only", "dada2_only", "both"}:
        return False, f"Unsupported analysis mode '{analysis_mode}'."
    return True, ""
