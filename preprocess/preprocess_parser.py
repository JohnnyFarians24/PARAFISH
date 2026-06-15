"""Parsers for preprocess outputs (trimmed and merged FASTQs)."""

from __future__ import annotations

import os
from pathlib import Path


def parse_preprocess_outputs(output_dir: str | Path) -> dict:
    """Build a normalized summary of preprocess output artifacts."""
    output_dir = Path(output_dir)
    if not output_dir.exists():
        return {"output_dir": str(output_dir), "samples": {}, "error": "Directory not found"}

    summary = {"output_dir": str(output_dir), "samples": {}}
    for sample_dir in sorted(output_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        collected = collect_preprocess_outputs(sample_dir)
        if collected:
            summary["samples"][sample_dir.name] = collected
    return summary


def collect_preprocess_outputs(output_dir: str | Path) -> dict:
    """Return available trimmed and merged files from a preprocess output directory."""
    output_dir = Path(output_dir)
    if not output_dir.exists():
        return {}

    outputs: dict[str, list[str]] = {"trimmed": [], "merged": []}
    for root, _, files in os.walk(output_dir):
        for name in files:
            lower = name.lower()
            path = os.path.join(root, name)
            if ".assembled.fastq" in lower:
                outputs["merged"].append(path)
            elif any(pattern in lower for pattern in ("_r1_trimmed", "_r2_trimmed", "_trimmed")) and lower.endswith((".fastq", ".fastq.gz")):
                outputs["trimmed"].append(path)

    return outputs if (outputs["trimmed"] or outputs["merged"]) else {}
