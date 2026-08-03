"""PEAR execution helpers for paired-end merging."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from core.config import PEAR_PATH
from core.models import AnalysisBranch, RunConfig, Sample, StepOutput


def _count_fastq_reads(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    line_count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for _ in handle:
            line_count += 1
    return line_count // 4


def run_pear(sample: Sample, config: RunConfig, trimmed_output: StepOutput) -> StepOutput:
    if config.branch != AnalysisBranch.VSEARCH:
        raise RuntimeError("PEAR is only used on the VSEARCH branch. DADA2 merges reads internally.")

    r1_trimmed = trimmed_output.files.get("r1_trimmed")
    r2_trimmed = trimmed_output.files.get("r2_trimmed")
    if r1_trimmed is None or r2_trimmed is None:
        return StepOutput(
            step_name="pear",
            success=False,
            files={},
            error_message="PEAR requer r1_trimmed e r2_trimmed no StepOutput do Trimmomatic.",
        )

    output_dir = config.output_dir / "preprocess" / sample.internal_id / "pear"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / sample.internal_id

    assembled = Path(f"{output_prefix}.assembled.fastq")
    # PEAR standard output is uncompressed; some versions may emit .assembled.fastq.gz.
    unassembled_fwd = Path(f"{output_prefix}.unassembled.forward.fastq")
    unassembled_rev = Path(f"{output_prefix}.unassembled.reverse.fastq")
    discarded = Path(f"{output_prefix}.discarded.fastq")
    log_path = output_dir / "pear.log"

    params = config.pear_params
    command = [
        PEAR_PATH,
        "-f",
        str(r1_trimmed),
        "-r",
        str(r2_trimmed),
        "-o",
        str(output_prefix),
        "-v",
        str(params.get("min_overlap", 10)),
        "-n",
        str(params.get("min_length", 200)),
        "-m",
        str(params.get("max_length", 500)),
        "-p",
        str(params.get("p_value", 0.01)),
        "-u",
        str(params.get("max_uncalled_base", 0.0)),
        "-j",
        str(params.get("threads", 4)),
        "-g",
        str(params.get("test", 1)),
    ]

    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        return StepOutput(
            step_name="pear",
            success=False,
            files={
                "assembled": assembled,
                "unassembled_fwd": unassembled_fwd,
                "unassembled_rev": unassembled_rev,
                "discarded": discarded,
                "log": log_path,
            },
            log=log_path,
            error_message=f"PEAR failed: {exc}",
        )

    log_path.write_text((completed.stdout or "") + "\n" + (completed.stderr or ""), encoding="utf-8")
    log_text = log_path.read_text(encoding="utf-8", errors="replace")

    assembled_reads = _count_fastq_reads(assembled)
    if assembled_reads == 0:
        return StepOutput(
            step_name="pear",
            success=False,
            files={
                "assembled": assembled,
                "unassembled_fwd": unassembled_fwd,
                "unassembled_rev": unassembled_rev,
                "discarded": discarded,
                "log": log_path,
            },
            log=log_path,
            error_message="PEAR finished but assembled.fastq is empty or missing.",
        )

    merge_rate = 0.0
    total_reads = 0.0
    discarded_rate = 0.0
    match = re.search(r"Assembled reads\s*[:=\s]+([0-9.]+)%", log_text, flags=re.IGNORECASE)
    if match:
        merge_rate = float(match.group(1)) / 100.0
    match = re.search(r"Total reads\s*[:=\s]+([0-9.]+)", log_text, flags=re.IGNORECASE)
    if match:
        total_reads = float(match.group(1))
    match = re.search(r"Discarded reads\s*[:=\s]+([0-9.]+)%", log_text, flags=re.IGNORECASE)
    if match:
        discarded_rate = float(match.group(1)) / 100.0

    warning = None
    if merge_rate < 0.50:
        warning = "PEAR merge rate below 50%."

    return StepOutput(
        step_name="pear",
        success=True,
        files={
            "assembled": assembled,
            "unassembled_fwd": unassembled_fwd,
            "unassembled_rev": unassembled_rev,
            "discarded": discarded,
            "log": log_path,
        },
        log=log_path,
        metrics={
            "merge_rate": merge_rate,
            "discarded_rate": discarded_rate,
            "total_reads": total_reads,
            "assembled_reads": float(assembled_reads),
        },
        error_message=warning,
    )

