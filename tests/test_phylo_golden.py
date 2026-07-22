"""The phylo golden contract: per-set trees the reader resolves, frozen byte-for-byte.

Phylo is the last stage into the golden net. Its trees are byte-reproducible — MAFFT and
FastTree are deterministic and RAxML-NG is pinned by ``seed=42`` — so the alignment and
both trees are frozen; the ``*.log``/``*.rba`` churn (wall-clock, binary caches) is not,
the same non-content the golden already excludes for FastQC zips and table.biom's date.

Two proofs, anchored on the COMMITTED golden so they run on a clean clone (not runs/,
which is gitignored):
 * the reader resolves the per-set layout the writer produces — the tree lives under its
   runner's subdir (``raxml/``, ``fasttree/``), NEVER flat at ``phylo/<set>/tree.nwk``,
   which is exactly where the pre-fix readers looked and found nothing ("árvore não
   encontrada");
 * the frozen support tree is a valid Newick whose tip count matches the alignment.

OTU-only by nature: fixture_asv has 3 ASVs and the aligner floor is 4 tips, so the ASV set
cannot exercise phylo at all. The golden pins the OTU/HYBRID path and leaves the ASV set
and the FAST/ROBUST-only modes UNexercised — the phylo analogue of the ASV-vs-aligner gap.
"""

from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# golden/<dataset>/<stage>/<rel>: the phylo stage freezes rels that start with "phylo/",
# so the frozen subtree is golden/fixture/phylo/phylo/otu/...
GOLDEN_PHYLO = ROOT / "golden" / "fixture" / "phylo" / "phylo"


class PhyloReaderResolutionTests(unittest.TestCase):
    """The reader finds the per-set tree the writer wrote — under its runner subdir."""

    def setUp(self):
        if not (GOLDEN_PHYLO / "otu").exists():
            self.skipTest("phylo golden not frozen")
        self._tmp = Path(tempfile.mkdtemp())
        shutil.copytree(GOLDEN_PHYLO, self._tmp / "phylo")

    def tearDown(self):
        shutil.rmtree(self._tmp, ignore_errors=True)

    def test_the_reader_resolves_the_per_set_support_tree(self):
        from ui.feature_sets import phylo_tree_path

        tree = phylo_tree_path(self._tmp, "otu")
        self.assertIsNotNone(tree)
        self.assertTrue(Path(tree).exists())
        # The bootstrap support tree wins — it is what the viewer colours by.
        self.assertTrue(tree.endswith("raxml.raxml.support"))

    def test_the_flat_paths_the_old_readers_used_do_not_exist(self):
        # Regression guard for "árvore não encontrada": the pre-fix readers looked at these
        # flat names, neither of which the writer ever produced.
        self.assertFalse((self._tmp / "phylo" / "otu" / "tree.nwk").exists())
        self.assertFalse((self._tmp / "phylo" / "tree.nwk").exists())

    def test_alignment_and_tip_mapping_resolve_per_set(self):
        from ui.feature_sets import phylo_artifact

        aligned = phylo_artifact(self._tmp, "otu", "aligned_fasta")
        mapping = phylo_artifact(self._tmp, "otu", "reversed_asvs_csv")
        self.assertTrue(aligned and Path(aligned).exists())
        self.assertTrue(mapping and Path(mapping).exists())

    def test_the_asv_set_has_no_tree(self):
        # fixture_asv is 3 tips < the 4-tip aligner floor, so nothing to resolve.
        from ui.feature_sets import phylo_tree_path

        self.assertIsNone(phylo_tree_path(self._tmp, "asv"))


class PhyloGoldenConsistencyTests(unittest.TestCase):
    """The frozen support tree is a valid Newick whose tip count matches the alignment."""

    def test_support_tree_tip_count_matches_the_alignment(self):
        if not (GOLDEN_PHYLO / "otu").exists():
            self.skipTest("phylo golden not frozen")
        support = (GOLDEN_PHYLO / "otu" / "raxml" / "raxml.raxml.support").read_text()
        aligned = (GOLDEN_PHYLO / "otu" / "mafft" / "aligned.fasta").read_text()
        # A Newick with N tips carries N-1 commas at top level; the alignment has N records.
        self.assertEqual(support.count(",") + 1, aligned.count(">"))


if __name__ == "__main__":
    unittest.main()
