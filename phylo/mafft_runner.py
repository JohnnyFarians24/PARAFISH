"""MAFFT execution helpers for sequence alignment."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Dict

from core.config import MAFFT_PATH
from core.models import StepOutput
from core.validation import validate_fasta_headers
from phylo.fasta_utils import sanitize_fasta_for_phylo


def run_mafft(fasta_path: Path, output_dir: Path, params: Dict) -> StepOutput:
    output_dir.mkdir(parents=True, exist_ok=True)
    sanitized = output_dir / "sanitized.fasta"
    # Name the tips after what they are: VSEARCH clusters are OTUs, DADA2 variants
    # are ASVs. ``feature_prefix`` is set by the phylo runner from the analysis method.
    mapping = sanitize_fasta_for_phylo(fasta_path, sanitized, prefix=str(params.get("feature_prefix") or "ASV"))

    # validate_fasta_headers returns a list of bad headers (truthy on errors)
    if validate_fasta_headers(sanitized):
        return StepOutput(
            step_name="mafft",
            success=False,
            files={},
            error_message="Invalid headers after sanitisation.",
        )

    aligned_fasta = output_dir / "aligned.fasta"
    log_path = output_dir / "mafft.log"
    reversed_csv = output_dir / "reversed_asvs.csv"

    threads = str(params.get("threads", 8))
    adjust = params.get("adjust_direction", True)

    primary = [MAFFT_PATH]
    if params.get("mode", "auto") == "accurate":
        primary += ["--maxiterate", "1000"]
    else:
        primary += ["--auto"]
    if adjust:
        primary.append("--adjustdirection")
    primary += ["--thread", threads, str(sanitized)]

    # Memory-frugal fallback: a single progressive pass without FFT uses far less
    # RAM, so alignment still completes on small machines where --auto is
    # OOM-killed (MAFFT then exits non-zero).
    fallback = [MAFFT_PATH, "--retree", "1", "--nofft"]
    if adjust:
        fallback.append("--adjustdirection")
    fallback += ["--thread", threads, str(sanitized)]

    def _attempt(command):
        return subprocess.run(command, capture_output=True, text=True, check=True)

    try:
        completed = _attempt(primary)
        aligned_fasta.write_text(completed.stdout or "", encoding="utf-8")
        log_path.write_text((completed.stderr or "") + "\n" + (completed.stdout or ""), encoding="utf-8")
    except subprocess.CalledProcessError as primary_exc:
        try:
            completed = _attempt(fallback)
            aligned_fasta.write_text(completed.stdout or "", encoding="utf-8")
            log_path.write_text(
                "MAFFT --auto failed (most likely out of memory); fell back to the "
                "economica --retree 1 --nofft.\n\n" + (completed.stderr or "") + "\n" + (completed.stdout or ""),
                encoding="utf-8",
            )
        except subprocess.CalledProcessError as fallback_exc:
            log_path.write_text(
                f"MAFFT failed with both strategies.\nPrimary: {primary_exc}\nMemory-saving: {fallback_exc}\n"
                "Likely cause: not enough memory (RAM) for this number of sequences.",
                encoding="utf-8",
            )
            return StepOutput(
                step_name="mafft",
                success=False,
                files={},
                log=log_path,
                error_message="MAFFT failed (most likely out of memory). See phylo/mafft/mafft.log.",
            )

    # write mapping of ASV -> original header
    with reversed_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ASV_ID", "Original_Header"])
        for asv_id, original in mapping.items():
            writer.writerow([asv_id, original])

    return StepOutput(
        step_name="mafft",
        success=True,
        files={
            "aligned_fasta": aligned_fasta,
            "log": log_path,
            "reversed_asvs_csv": reversed_csv,
        },
        log=log_path,
    )
