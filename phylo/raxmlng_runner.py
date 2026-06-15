"""RAxML-NG execution helpers for phylogenetic inference."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict

from core.config import RAXML_PATH
from core.models import StepOutput
from core.validation import validate_fasta_headers


def run_raxml(aligned_fasta: Path, output_dir: Path, params: Dict) -> StepOutput:
    output_dir.mkdir(parents=True, exist_ok=True)

    errors = validate_fasta_headers(aligned_fasta)
    if errors:
        return StepOutput(
            step_name="raxml",
            success=False,
            files={},
            error_message=f"Invalid FASTA headers: {'; '.join(errors[:3])}",
        )

    prefix = output_dir / "raxml"
    log_path = output_dir / "raxml.log"
    model = params.get("model", "GTR+G")
    threads = params.get("threads", 8)
    seed = params.get("seed", 42)
    bs_metric = params.get("bs_metric", "fbp")
    outgroup = params.get("outgroup")
    search_tree = params.get("search_tree")

    commands = [
        [
            RAXML_PATH, "--parse",
            "--msa", str(aligned_fasta),
            "--model", str(model),
            "--prefix", f"{prefix}_PARSE",
        ],
        [
            RAXML_PATH, "--all",
            "--msa", str(aligned_fasta),
            "--model", str(model),
            "--prefix", str(prefix),
            "--threads", f"auto{{{threads}}}",
            "--seed", str(seed),
            "--bs-metric", str(bs_metric),
            "--bs-trees", str(int(params.get("bs_trees", 10))),
            "--redo",
        ]
        + (["--outgroup", str(outgroup)] if outgroup else [])
            + (["--tree", str(search_tree)] if search_tree and Path(search_tree).exists() else []),
    ]

    log_lines: list[str] = []
    try:
        for command in commands:
            completed = subprocess.run(command, capture_output=True, text=True, check=True)
            log_lines.append(completed.stdout or "")
            log_lines.append(completed.stderr or "")
    except subprocess.CalledProcessError as exc:
        log_lines.append(str(exc))
        log_path.write_text("\n".join(log_lines), encoding="utf-8")
        return StepOutput(
            step_name="raxml",
            success=False,
            files={},
            log=log_path,
            error_message=f"RAxML-NG falhou: {exc}",
        )

    log_path.write_text("\n".join(log_lines), encoding="utf-8")

    bs_suffix = "supportFBP" if bs_metric == "fbp" else "supportTBE"
    best_tree = output_dir / "raxml.raxml.bestTree"
    support_tree = output_dir / f"raxml.raxml.{bs_suffix}"
    bootstraps = output_dir / "raxml.raxml.bootstraps"
    best_model = output_dir / "raxml.raxml.bestModel"

    if not best_tree.exists():
        return StepOutput(
            step_name="raxml",
            success=False,
            files={"log": log_path},
            log=log_path,
            error_message="RAxML-NG terminou mas bestTree não foi gerado. Verifica o log.",
        )

    return StepOutput(
        step_name="raxml",
        success=True,
        files={
            "best_tree": best_tree,
            "support_tree": support_tree,
            "bootstraps": bootstraps,
            "best_model": best_model,
            "log": log_path,
        },
        log=log_path,
    )
