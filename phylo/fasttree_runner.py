"""FastTree execution helpers for tree inference."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict

from core.config import FASTTREE_PATH
from core.models import StepOutput
from core.validation import validate_fasta_headers
from phylo import MIN_TREE_SEQUENCES


def run_fasttree(aligned_fasta: Path, output_dir: Path, params: Dict) -> StepOutput:
    output_dir.mkdir(parents=True, exist_ok=True)

    if validate_fasta_headers(aligned_fasta):
        return StepOutput(
            step_name="fasttree",
            success=False,
            files={},
            error_message="Alignment contains invalid headers.",
        )

    sequence_count = sum(
        1 for line in aligned_fasta.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(">")
    )
    if sequence_count < MIN_TREE_SEQUENCES:
        return StepOutput(
            step_name="fasttree",
            success=False,
            files={},
            error_message=f"Alignment has only {sequence_count} sequences (minimum: {MIN_TREE_SEQUENCES}).",
        )

    tree_nwk = output_dir / "tree.nwk"
    log_path = output_dir / "fasttree.log"

    command = [FASTTREE_PATH, "-nt"]
    if params.get("gtr", True):
        command.append("-gtr")
    if params.get("gamma", True):
        command.append("-gamma")
    if params.get("fastest", False):
        command.append("-fastest")
    if params.get("no_support", False):
        command.append("-nosupport")
    command.append(str(aligned_fasta))

    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=True)
        tree_nwk.write_text(completed.stdout or "", encoding="utf-8")
        log_path.write_text(completed.stderr or "", encoding="utf-8")
    except subprocess.CalledProcessError as exc:
        log_path.write_text(str(exc), encoding="utf-8")
        return StepOutput(
            step_name="fasttree",
            success=False,
            files={},
            log=log_path,
            error_message=f"FastTree failed: {exc}",
        )

    if not tree_nwk.exists() or tree_nwk.stat().st_size == 0:
        return StepOutput(
            step_name="fasttree",
            success=False,
            files={"log": log_path},
            log=log_path,
            error_message="FastTree finished without producing a tree. Check the input alignment.",
        )

    return StepOutput(
        step_name="fasttree",
        success=True,
        files={"tree_nwk": tree_nwk, "log": log_path},
        log=log_path,
    )
