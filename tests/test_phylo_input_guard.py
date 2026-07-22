"""The phylo input guard must refuse at the SAME floor the aligner enforces.

The honesty bug: the UI told the user two sequences were enough for a tree (``seq_count < 2``),
but FastTree — which the default HYBRID and FAST modes both route through — refuses below four
("Alinhamento com apenas N sequências (mínimo: 4)"). So a 3-tip run passed the UI guard, waited
through QC + preprocess + analysis, and only failed DOWNSTREAM at the aligner — a false "you may
proceed" of the "server did not respond" family. The floor of 4 is measured, not guessed: the
fixture_asv set has 3 ASVs and hits it.

The fix is at the source of the divergence: ONE floor (``MIN_TREE_SEQUENCES``), shared by the
tree runner and both UI guards (the manual button and the pipeline), so the number the user is
refused at cannot drift from the number the aligner enforces. These tests lock that a 3-tip input
is refused AT THE ENTRANCE, named with the real floor.
"""

from __future__ import annotations

import inspect
import tempfile
import unittest
from pathlib import Path


class PhyloEntryGuardTests(unittest.TestCase):
    def test_the_entry_guard_refuses_three_tips_and_names_the_real_floor(self):
        from phylo import MIN_TREE_SEQUENCES, phylo_input_error

        message = phylo_input_error(3)
        self.assertIsNotNone(message, "3 tips must be refused at the entry guard, not downstream")
        self.assertIn(str(MIN_TREE_SEQUENCES), message)
        self.assertIn("4", message)
        self.assertNotIn("mínimo 2", message)

    def test_a_sufficient_input_passes_the_entry_guard(self):
        from phylo import MIN_TREE_SEQUENCES, phylo_input_error

        self.assertIsNone(phylo_input_error(MIN_TREE_SEQUENCES))
        self.assertIsNone(phylo_input_error(10))

    def test_the_guard_floor_equals_the_aligner_floor(self):
        # The whole bug was the guard (2) and the aligner (4) disagreeing. Pin them equal:
        # the same 3-tip input the guard now refuses is exactly what FastTree refuses.
        from phylo import MIN_TREE_SEQUENCES
        from phylo.fasttree_runner import run_fasttree

        self.assertEqual(MIN_TREE_SEQUENCES, 4)
        with tempfile.TemporaryDirectory() as d:
            aligned = Path(d) / "aligned.fasta"
            aligned.write_text(">A\nACGT\n>B\nACGT\n>C\nACGT\n", encoding="utf-8")
            result = run_fasttree(aligned, Path(d) / "out", {})
            self.assertFalse(result.success)
            self.assertIn(str(MIN_TREE_SEQUENCES), result.error_message)

    def test_both_ui_guards_use_the_shared_floor(self):
        # Guard against the helper existing but a guard site still hard-coding its own number:
        # both the manual button and the pipeline must resolve the floor from the shared source.
        import ui.callbacks_phylo as manual
        import ui.callbacks_progress as pipeline

        for module in (manual, pipeline):
            source = inspect.getsource(module)
            self.assertNotIn("mínimo 2", source, f"{module.__name__} still hard-codes the old floor of 2")
            self.assertTrue(
                "phylo_input_error" in source or "MIN_TREE_SEQUENCES" in source,
                f"{module.__name__} does not resolve the phylo floor from the shared source",
            )


if __name__ == "__main__":
    unittest.main()
