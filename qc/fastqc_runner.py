"""FastQC command orchestration helpers."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path

from core.config import FASTQC_PATH
from core.models import AnalysisMode, Sample, StepOutput


def _fastqc_base_name(input_file: Path) -> str:
    name = input_file.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return input_file.stem


def _fastqc_report_paths(input_file: str | Path, output_dir: str | Path):
    input_path = Path(input_file)
    output_path = Path(output_dir)
    stem = _fastqc_base_name(input_path)
    return {
        "html": output_path / f"{stem}_fastqc.html",
        "zip": output_path / f"{stem}_fastqc.zip",
    }


def _parse_fastqc_zip(zip_path: Path, suffix: str) -> dict[str, float]:
    metrics: dict[str, float] = {}
    try:
        with zipfile.ZipFile(zip_path) as archive:
            data_files = [name for name in archive.namelist() if name.endswith("fastqc_data.txt")]
            if not data_files:
                return metrics
            content = archive.read(data_files[0]).decode("utf-8", errors="replace")
        for line in content.splitlines():
            if line.startswith("Total Sequences"):
                parts = line.split("\t", 1)
                if len(parts) == 2:
                    metrics[f"total_sequences_{suffix}"] = float(parts[1])
            elif line.startswith("%GC"):
                parts = line.split("\t", 1)
                if len(parts) == 2:
                    metrics[f"pct_gc_{suffix}"] = float(parts[1])
            elif line and not line.startswith("#") and "\tFAIL" in line:
                key = f"failed_modules_{suffix}"
                metrics[key] = metrics.get(key, 0.0) + 1.0
    except (zipfile.BadZipFile, KeyError, UnicodeDecodeError, OSError, ValueError):
        return metrics
    metrics.setdefault(f"failed_modules_{suffix}", 0.0)
    return metrics


def run_fastqc(sample: Sample, output_dir: Path, threads: int = 4) -> StepOutput:
    """Run FastQC on the raw reads for a single sample and return a StepOutput."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / f"{sample.internal_id}_fastqc.log"
    inputs = [sample.paths.r1]
    if sample.mode == AnalysisMode.PE and sample.paths.r2:
        inputs.append(sample.paths.r2)
    fastqc_threads = 1 if (sample.mode == AnalysisMode.PE and sample.paths.r2) else min(threads, 4)

    command = [FASTQC_PATH]
    command.extend(str(path) for path in inputs)
    command.extend(["--outdir", str(output_dir), "--threads", str(fastqc_threads), "--quiet"])

    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=True)
        log_path.write_text((completed.stdout or "") + ("\n" if completed.stdout and completed.stderr else "") + (completed.stderr or ""), encoding="utf-8")
    except FileNotFoundError as exc:
        log_path.write_text(str(exc), encoding="utf-8")
        return StepOutput(
            step_name="fastqc",
            success=False,
            files={},
            log=log_path,
            error_message=f"FastQC not found: {exc}",
        )
    except subprocess.CalledProcessError as exc:
        log_path.write_text((exc.stdout or "") + ("\n" if exc.stdout and exc.stderr else "") + (exc.stderr or ""), encoding="utf-8")
        return StepOutput(
            step_name="fastqc",
            success=False,
            files={},
            log=log_path,
            error_message=f"FastQC failed: {str(exc.stderr or exc.stdout or exc)[:500]}",
        )

    files: dict[str, Path] = {}
    metrics: dict[str, float] = {}

    r1_reports = _fastqc_report_paths(sample.paths.r1, output_dir)
    if r1_reports["html"].exists():
        files["html_r1"] = r1_reports["html"]
    if r1_reports["zip"].exists():
        files["zip_r1"] = r1_reports["zip"]
        metrics.update(_parse_fastqc_zip(r1_reports["zip"], "r1"))

    if sample.mode == AnalysisMode.PE and sample.paths.r2:
        r2_reports = _fastqc_report_paths(sample.paths.r2, output_dir)
        if r2_reports["html"].exists():
            files["html_r2"] = r2_reports["html"]
        if r2_reports["zip"].exists():
            files["zip_r2"] = r2_reports["zip"]
            metrics.update(_parse_fastqc_zip(r2_reports["zip"], "r2"))

    if not files:
        return StepOutput(
            step_name="fastqc",
            success=False,
            files={},
            log=log_path,
            error_message="FastQC exited cleanly but produced none of the expected reports. Check the outdir and the input file names.",
        )

    return StepOutput(
        step_name="fastqc",
        success=True,
        files=files,
        log=log_path,
        metrics=metrics,
    )
