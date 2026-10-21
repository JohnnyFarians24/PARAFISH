"""Validation helpers for contracts and filesystem-backed pipeline inputs."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import List

from .models import AnalysisBranch, AnalysisMode, RunConfig, Sample


INVALID_SAMPLE_NAME_CHARS = {" ", ":", ",", "(", ")", "/"}
INVALID_FASTA_HEADER_CHARS = {":", ",", "(", ")", " "}


def _count_reads_with_wc(path: Path) -> int:
    if str(path).endswith(".gz"):
        return -1
    result = subprocess.run(["wc", "-l", str(path)], capture_output=True, text=True, check=True)
    return int(result.stdout.strip().split()[0])


def validate_sample_names(samples: List[Sample]) -> List[str]:
    """Ensure no internal_id contains a forbidden character."""
    errors: List[str] = []
    for sample in samples:
        for char in INVALID_SAMPLE_NAME_CHARS:
            if char in sample.internal_id:
                errors.append(f"Invalid sample ID '{sample.internal_id}': contains '{char}'.")
                break
    return errors


def validate_manifest(config: RunConfig) -> List[str]:
    """Validate sample layout, duplicate IDs and branch-specific requirements."""
    errors: List[str] = []
    seen_ids = set()
    for sample in config.samples:
        if sample.internal_id in seen_ids:
            errors.append(f"Duplicate sample ID: {sample.internal_id}")
        seen_ids.add(sample.internal_id)

        try:
            sample.paths.validate(sample.mode)
        except AssertionError as exc:
            errors.append(str(exc))
            continue

        if sample.mode == AnalysisMode.PE:
            if sample.paths.r2 is None:
                errors.append(f"Modo PE requer R2: {sample.internal_id}")
            else:
                try:
                    reads_r1 = _count_reads_with_wc(sample.paths.r1)
                    reads_r2 = _count_reads_with_wc(sample.paths.r2)
                    if reads_r1 != -1 and reads_r2 != -1 and reads_r1 != reads_r2:
                        errors.append(f"R1 and R2 have a different number of lines for {sample.internal_id}.")
                except Exception as exc:
                    errors.append(f"Failed to count reads for {sample.internal_id}: {exc}")

        if sample.mode == AnalysisMode.SE and not sample.paths.r1.exists():
            errors.append(f"Modo SE requer R1: {sample.internal_id}")

    if getattr(config, "analysis_method", "vsearch") in {"dada2", "both"} and config.dada2_params.get("ref_fasta") is None:
        errors.append("DADA2 requires a reference database (ref_fasta). Set the path in dada2_params.")

    return errors


def validate_fasta_headers(fasta_path: Path) -> List[str]:
    """Verifies that no header contains prohibited characters."""
    errors: List[str] = []
    if not fasta_path.exists():
        return [f"FASTA file not found: {fasta_path}"]
    with fasta_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if any(char in header for char in INVALID_FASTA_HEADER_CHARS):
                errors.append(header)
    return errors
