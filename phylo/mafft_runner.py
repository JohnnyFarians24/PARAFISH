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
    mapping = sanitize_fasta_for_phylo(fasta_path, sanitized)

    # validate_fasta_headers returns a list of bad headers (truthy on errors)
    if validate_fasta_headers(sanitized):
        return StepOutput(
            step_name="mafft",
            success=False,
            files={},
            error_message="Cabeçalhos inválidos após sanitização.",
        )

    aligned_fasta = output_dir / "aligned.fasta"
    log_path = output_dir / "mafft.log"
    reversed_csv = output_dir / "reversed_asvs.csv"

    command = [MAFFT_PATH]
    if params.get("mode", "auto") == "accurate":
        command += ["--maxiterate", "1000"]
    else:
        command += ["--auto"]
    if params.get("adjust_direction", True):
        command.append("--adjustdirection")
    command += ["--thread", str(params.get("threads", 8)), str(sanitized)]

    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=True)
        aligned_fasta.write_text(completed.stdout or "", encoding="utf-8")
        log_path.write_text(
            (completed.stderr or "") + "\n" + (completed.stdout or ""),
            encoding="utf-8",
        )
    except subprocess.CalledProcessError as exc:
        log_path.write_text(str(exc), encoding="utf-8")
        return StepOutput(
            step_name="mafft",
            success=False,
            files={},
            log=log_path,
            error_message=f"MAFFT falhou: {exc}",
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
