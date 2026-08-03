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


def _finalize_phylo_output(result: StepOutput) -> StepOutput:
    """Expose a canonical ``tree_nwk`` key regardless of the runner used.

    FastTree emits ``tree_nwk`` while RAxML-NG emits ``best_tree`` /
    ``support_tree``. Downstream viewers (Results preview, phylo panel) expect a
    single primary Newick key, so alias one here. For RAxML-NG we prefer the
    bootstrap ``support_tree`` (it carries branch-support values) over the plain
    ``best_tree`` so the surfaced/visualized tree shows bootstrap support.
    """
    if result.success and isinstance(result.files, dict) and "tree_nwk" not in result.files:
        support = result.files.get("support_tree")
        best = result.files.get("best_tree")
        # Prefer the support tree only when it was actually produced; otherwise
        # fall back to the best tree so a missing bootstrap file never hides the
        # tree from downstream viewers.
        primary = support if support is not None and Path(support).exists() else best
        if primary is not None:
            result.files["tree_nwk"] = primary
    return result


def _with_alignment_files(result: StepOutput, mafft_out: StepOutput) -> StepOutput:
    """Carry MAFFT's artifacts onto the final output.

    ``run_phylo`` used to return only the tree builder's StepOutput, so the
    alignment and the ASV-id mapping never reached the UI even though both sit on
    disk — the alignment QC panel and Faith's PD id mapping need them.
    """
    if isinstance(result.files, dict):
        for key in ("aligned_fasta", "reversed_asvs_csv"):
            value = mafft_out.files.get(key)
            if value is not None:
                result.files.setdefault(key, value)
    return result


def run_phylo(aligned_input_fasta: Path, config: RunConfig, output_dir: Path) -> StepOutput:
    """Orchestrate the phylogenetic pipeline according to the configured mode."""
    output_dir.mkdir(parents=True, exist_ok=True)

    aligned_input = Path(aligned_input_fasta)
    if not aligned_input.exists():
        return StepOutput(
            step_name="phylo",
            success=False,
            files={},
            error_message=f"Input alignment not found: {aligned_input}.",
        )

    # VSEARCH clusters are OTUs; DADA2 variants are ASVs. Labelling clusters "ASV"
    # misdescribes what the tree's tips are.
    mafft_params = dict(config.mafft_params)
    mafft_params.setdefault(
        "feature_prefix",
        "ASV" if str(getattr(config, "analysis_method", "vsearch")).lower() == "dada2" else "OTU",
    )
    mafft_out = run_mafft(aligned_input, output_dir / "mafft", mafft_params)
    if not mafft_out.success:
        return mafft_out

    aligned = mafft_out.files["aligned_fasta"]

    if config.phylo_mode == PhyloMode.FAST:
        return _with_alignment_files(_finalize_phylo_output(run_fasttree(aligned, output_dir / "fasttree", config.phylo_params)), mafft_out)

    if config.phylo_mode == PhyloMode.ROBUST:
        return _with_alignment_files(_finalize_phylo_output(run_raxml(aligned, output_dir / "raxml", config.phylo_params)), mafft_out)

    # PhyloMode.AUTO: FastTree first, then RAxML seeded with that tree
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
            error_message="FastTree returned no tree_nwk to seed RAxML-NG with.",
        )
    if not tree_seed.exists():
        return StepOutput(
            step_name="phylo",
            success=False,
            files={"log": fasttree_out.files.get("log") or output_dir / "fasttree" / "fasttree.log"},
            error_message=f"FastTree devolveu um seed inexistente: {tree_seed}.",
        )
    raxml_params.setdefault("search_tree", str(tree_seed))
    return _with_alignment_files(_finalize_phylo_output(run_raxml(aligned, output_dir / "raxml", raxml_params)), mafft_out)
