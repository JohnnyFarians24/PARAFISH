"""Phylogenetic pipeline orchestrator."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict

from core.models import PhyloMode, RunConfig, StepOutput
from phylo.fasta_utils import sanitize_fasta_for_phylo
from phylo.mafft_runner import run_mafft
from phylo.fasttree_runner import run_fasttree
from phylo.raxmlng_runner import run_raxml


def run_phylo(aligned_input_fasta: Path, config: RunConfig, output_dir: Path) -> StepOutput:
    """Orchestrate the phylogenetic pipeline according to the configured mode."""
    output_dir.mkdir(parents=True, exist_ok=True)

    aligned_input = Path(aligned_input_fasta)
    if not aligned_input.exists():
        return StepOutput(
            step_name="phylo",
            success=False,
            files={},
            error_message=f"Alinhamento de entrada não encontrado: {aligned_input}.",
        )

    mafft_out = run_mafft(aligned_input, output_dir / "mafft", config.mafft_params)
    if not mafft_out.success:
        return mafft_out

    aligned = mafft_out.files["aligned_fasta"]

    if config.phylo_mode == PhyloMode.FAST:
        return run_fasttree(aligned, output_dir / "fasttree", config.phylo_params)

    if config.phylo_mode == PhyloMode.ROBUST:
        return run_raxml(aligned, output_dir / "raxml", config.phylo_params)

    # PhyloMode.AUTO: FastTree primeiro, depois RAxML com a árvore como seed
    fasttree_out = run_fasttree(aligned, output_dir / "fasttree", config.phylo_params)
    if not fasttree_out.success:
        return fasttree_out

    raxml_params = dict(config.phylo_params)
    tree_seed = fasttree_out.files.get("tree_nwk")
    if not isinstance(tree_seed, Path):
        return StepOutput(
            step_name="phylo",
            success=False,
            files={"log": fasttree_out.files.get("log") or output_dir / "fasttree" / "fasttree.log"},
            error_message="FastTree não devolveu tree_nwk para usar como seed do RAxML-NG.",
        )
    if not tree_seed.exists():
        return StepOutput(
            step_name="phylo",
            success=False,
            files={"log": fasttree_out.files.get("log") or output_dir / "fasttree" / "fasttree.log"},
            error_message=f"FastTree devolveu um seed inexistente: {tree_seed}.",
        )
    raxml_params.setdefault("search_tree", str(tree_seed))
    return run_raxml(aligned, output_dir / "raxml", raxml_params)
