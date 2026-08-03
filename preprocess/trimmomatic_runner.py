"""Trimmomatic execution helpers for preprocessing."""

from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from core.config import TRIMMOMATIC_ADAPTERS, TRIMMOMATIC_ADAPTERS_PE, TRIMMOMATIC_JAR_PATH
from core.models import AnalysisMode, RunConfig, Sample, StepOutput


def _count_fastq_reads(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    line_count = 0
    if path.suffix == ".gz":
        with path.open("rb") as probe:
            magic = probe.read(2)
        opener = gzip.open if magic == b"\x1f\x8b" else open
    else:
        opener = open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for _ in handle:
            line_count += 1
    return line_count // 4


def _build_output_prefix(sample: Sample, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / sample.internal_id


def _resolve_adapter_path(adapter: str) -> Path:
    candidate = Path(adapter)
    if candidate.exists():
        return candidate
    local_candidate = Path.cwd() / "Trimmomatic-0.39" / "adapters" / candidate.name
    return local_candidate


def _run_command(command: List[str], log_path: Path) -> subprocess.CompletedProcess:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    log_path.write_text((completed.stdout or "") + "\n" + (completed.stderr or ""), encoding="utf-8")
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, command, completed.stdout, completed.stderr)
    return completed


def run_trimmomatic(sample: Sample, config: RunConfig) -> StepOutput:
    """Run Trimmomatic for a single sample and return the standard step output."""
    output_dir = config.output_dir / "preprocess" / sample.internal_id
    prefix = _build_output_prefix(sample, output_dir)
    log_path = output_dir / "trimlog.txt"

    adapter = TRIMMOMATIC_ADAPTERS_PE if sample.mode == AnalysisMode.PE else TRIMMOMATIC_ADAPTERS
    adapter_path = _resolve_adapter_path(adapter)
    if not adapter_path.exists():
        return StepOutput(
            step_name="trimmomatic",
            success=False,
            files={},
            error_message=f"Adapter file not found: {adapter}. Set TRIMMOMATIC_ADAPTERS_PE in environment.",
        )

    if sample.mode == AnalysisMode.PE and (sample.paths.r2 is None or not sample.paths.r2.exists()):
        return StepOutput(
            step_name="trimmomatic",
            success=False,
            files={},
            log=None,
            error_message=f"Paired-end sample with no valid R2: {sample.internal_id}",
        )

    if sample.mode == AnalysisMode.PE:
        r1_out = Path(f"{prefix}_R1_trimmed.fastq.gz")
        r2_out = Path(f"{prefix}_R2_trimmed.fastq.gz")
        r1_unpaired = Path(f"{prefix}_R1_unpaired.fastq.gz")
        r2_unpaired = Path(f"{prefix}_R2_unpaired.fastq.gz")
        command = [
            "java",
            "-jar",
            TRIMMOMATIC_JAR_PATH,
            "PE",
            # Force Phred+33. Trimmomatic otherwise sniffs the encoding from the first
            # reads, and a file whose quality line is uniform (e.g. all '?') sits in the
            # ambiguous band between Phred+33 and +64, so detection fails with
            # "Unable to detect quality encoding" and the whole sample is lost. Every
            # Illumina run since ~2011 is Phred+33, so asserting it is safe and removes a
            # class of spurious failures that have nothing to do with the data's quality.
            "-phred33",
            str(sample.paths.r1),
            str(sample.paths.r2),
            str(r1_out),
            str(r1_unpaired),
            str(r2_out),
            str(r2_unpaired),
            f"ILLUMINACLIP:{adapter_path}:2:30:10",
            f"LEADING:{config.trimmomatic_params.get('leading', 3)}",
            f"TRAILING:{config.trimmomatic_params.get('trailing', 3)}",
            f"SLIDINGWINDOW:{config.trimmomatic_params.get('window_size', 4)}:{config.trimmomatic_params.get('window_quality', 20)}",
            f"MINLEN:{config.trimmomatic_params.get('minlen', 200)}",
        ]
        files = {"r1_trimmed": r1_out, "r2_trimmed": r2_out, "log": log_path}
    else:
        r1_out = Path(f"{prefix}_trimmed.fastq.gz")
        command = [
            "java",
            "-jar",
            TRIMMOMATIC_JAR_PATH,
            "SE",
            "-phred33",  # see the PE branch: avoids "Unable to detect quality encoding"
            str(sample.paths.r1),
            str(r1_out),
            f"ILLUMINACLIP:{adapter_path}:2:30:10",
            f"LEADING:{config.trimmomatic_params.get('leading', 3)}",
            f"TRAILING:{config.trimmomatic_params.get('trailing', 3)}",
            f"SLIDINGWINDOW:{config.trimmomatic_params.get('window_size', 4)}:{config.trimmomatic_params.get('window_quality', 20)}",
            f"MINLEN:{config.trimmomatic_params.get('minlen', 200)}",
        ]
        files = {"r1_trimmed": r1_out, "log": log_path}

    try:
        _run_command(command, log_path)
    except subprocess.CalledProcessError as exc:
        return StepOutput(
            step_name="trimmomatic",
            success=False,
            files=files,
            log=log_path,
            error_message=f"Trimmomatic failed: {exc}",
        )

    read_counts = []
    for key in ("r1_trimmed", "r2_trimmed"):
        file_path = files.get(key)
        if file_path and file_path.exists():
            read_counts.append(_count_fastq_reads(file_path))

    if not read_counts or max(read_counts) == 0:
        return StepOutput(
            step_name="trimmomatic",
            success=False,
            files=files,
            log=log_path,
            error_message="Trimming produced 0 reads. Check the MINLEN parameter and the adapter file.",
        )

    return StepOutput(
        step_name="trimmomatic",
        success=True,
        files=files,
        log=log_path,
        metrics={"trimmed_reads": float(max(read_counts))},
    )

