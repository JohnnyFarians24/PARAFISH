"""Phylogenetic processing module for alignment and tree inference."""

from __future__ import annotations

from typing import Optional

# The minimum number of sequences a tree can be built from. FastTree — which the default
# HYBRID and FAST modes both route through — refuses below four, because an unrooted tree
# needs at least four tips to be non-trivial. This is the ONE definition of that floor: the
# tree runner and the UI guards (manual button and pipeline) all read it, so the number the
# user is refused at cannot drift from the number the aligner enforces downstream. The bug
# this closes had the UI guard at 2 while the aligner was at 4, so a 3-tip run passed the
# guard and only failed after QC + preprocess + analysis had already run.
MIN_TREE_SEQUENCES = 4


def phylo_input_error(seq_count: int) -> Optional[str]:
    """The reason a FASTA is too small to build a tree, or ``None`` when it is big enough.

    Shared by both UI guards so the refusal happens AT THE ENTRANCE, named with the real
    floor, instead of blowing up downstream at the aligner.
    """
    if seq_count < MIN_TREE_SEQUENCES:
        return f"FASTA tem {seq_count} sequências; mínimo {MIN_TREE_SEQUENCES} para phylo"
    return None
