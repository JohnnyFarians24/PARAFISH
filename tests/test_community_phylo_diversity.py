"""Tests for community ecology + phylogenetic diversity numeric cores."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from analysis.community_analysis import (
    bray_curtis_matrix,
    collapsed_taxa_table,
    pcoa,
    sample_summary,
    taxa_composition,
)
from phylo.phylo_diversity import faith_pd, tree_layout


_MATRIX = "#OTU ID,S1,S2,S3\nOTU1,120,80,0\nOTU2,5,5,5\nOTU3,0,0,42\nOTU4,10,0,0\n"
_TAXONOMY = (
    "ASV_ID,Consensus_Kingdom,Consensus_Genus\n"
    "OTU1,Animalia,Gyrodactylus\n"
    "OTU2,Animalia,Dactylogyrus\n"
    "OTU3,Animalia,\n"
    "OTU4,Animalia,Gyrodactylus\n"
)


def _write(temp_path, text, name):
    p = temp_path / name
    p.write_text(text, encoding="utf-8")
    return p


class BetaDiversityTests(unittest.TestCase):
    def test_bray_curtis_known_values(self):
        with tempfile.TemporaryDirectory() as d:
            m = _write(Path(d), _MATRIX, "otu_matrix.csv")
            bc = bray_curtis_matrix(m)
            self.assertEqual(bc["samples"], ["S1", "S2", "S3"])
            mat = bc["matrix"]
            self.assertAlmostEqual(mat[0][1], 50 / 220, places=4)   # S1 vs S2
            self.assertAlmostEqual(mat[0][2], 172 / 182, places=4)  # S1 vs S3
            self.assertAlmostEqual(mat[1][2], 122 / 132, places=4)  # S2 vs S3
            self.assertEqual(mat[0][0], 0.0)

    def test_pcoa_separates_identical_from_different(self):
        # A and B identical; C and D identical but different community.
        matrix = "#OTU ID,A,B,C,D\nOTU1,100,100,0,0\nOTU2,0,0,100,100\n"
        with tempfile.TemporaryDirectory() as d:
            m = _write(Path(d), matrix, "otu_matrix.csv")
            result = pcoa(bray_curtis_matrix(m))
            coords = dict(zip(result["samples"], result["coords"]))
            # Identical samples land at (nearly) the same coordinate.
            self.assertAlmostEqual(coords["A"][0], coords["B"][0], places=4)
            self.assertAlmostEqual(coords["C"][0], coords["D"][0], places=4)
            # Different communities separate along PC1.
            self.assertNotAlmostEqual(coords["A"][0], coords["C"][0], places=3)

    def test_pcoa_degenerate_two_samples(self):
        matrix = "#OTU ID,A,B\nOTU1,10,0\n"
        with tempfile.TemporaryDirectory() as d:
            m = _write(Path(d), matrix, "otu_matrix.csv")
            result = pcoa(bray_curtis_matrix(m))
            self.assertEqual(len(result["coords"]), 2)


class CompositionTests(unittest.TestCase):
    def test_composition_collapses_by_genus(self):
        with tempfile.TemporaryDirectory() as d:
            m = _write(Path(d), _MATRIX, "otu_matrix.csv")
            t = _write(Path(d), _TAXONOMY, "taxonomy_consensus.csv")
            comp = taxa_composition(m, t, rank="genus", top_n=10)
            self.assertEqual(comp["samples"], ["S1", "S2", "S3"])
            gyro = comp["taxa"].index("Gyrodactylus")
            # Gyrodactylus = OTU1+OTU4: S1 = 130, relative 130/135.
            self.assertAlmostEqual(comp["rel"][gyro][0], 130 / 135, places=4)
            self.assertIn("Unclassified", comp["taxa"])  # OTU3 has no genus
            # Each sample's relative abundances sum to ~1.
            for j in range(3):
                self.assertAlmostEqual(sum(comp["rel"][t][j] for t in range(len(comp["taxa"]))), 1.0, places=4)

    def test_top_n_folds_into_other(self):
        rows = "\n".join(f"OTU{i},{100 - i},{100 - i}" for i in range(1, 11))
        matrix = "#OTU ID,S1,S2\n" + rows + "\n"
        tax = "ASV_ID,Consensus_Genus\n" + "\n".join(f"OTU{i},Genus{i}" for i in range(1, 11)) + "\n"
        with tempfile.TemporaryDirectory() as d:
            m = _write(Path(d), matrix, "otu_matrix.csv")
            t = _write(Path(d), tax, "taxonomy_consensus.csv")
            comp = taxa_composition(m, t, rank="genus", top_n=3)
            self.assertEqual(len(comp["taxa"]), 4)  # 3 + Other
            self.assertEqual(comp["taxa"][-1], "Other")

    def test_collapsed_table_written(self):
        with tempfile.TemporaryDirectory() as d:
            m = _write(Path(d), _MATRIX, "otu_matrix.csv")
            t = _write(Path(d), _TAXONOMY, "taxonomy_consensus.csv")
            out = collapsed_taxa_table(m, t, "genus", Path(d) / "taxa_table.csv")
            self.assertTrue(out and out.exists())
            content = out.read_text(encoding="utf-8")
            self.assertIn("Gyrodactylus", content)
            self.assertIn("S1", content)

    def test_sample_summary_pct_classified(self):
        with tempfile.TemporaryDirectory() as d:
            m = _write(Path(d), _MATRIX, "otu_matrix.csv")
            t = _write(Path(d), _TAXONOMY, "taxonomy_consensus.csv")
            summary = {s["sample"]: s for s in sample_summary(m, t, rank="genus")}
            # S1 reads=135, classified reads = OTU1(120)+OTU2(5)+OTU4(10)=135 (OTU3 absent in S1) -> 100%
            self.assertEqual(summary["S1"]["reads"], 135)
            self.assertEqual(summary["S1"]["pct_classified"], 100.0)
            # S3 reads=47, classified = OTU2(5); OTU3(42) unclassified -> ~10.6%
            self.assertAlmostEqual(summary["S3"]["pct_classified"], round(100 * 5 / 47, 1), places=1)


class FaithPDTests(unittest.TestCase):
    _TREE = "((A:1,B:1):1,(C:1,D:1):1);"
    _PD_MATRIX = "#OTU ID,S1,S2,S3\nA,5,5,5\nB,5,0,0\nC,0,5,0\nD,0,0,0\n"

    def test_faith_pd_rooted(self):
        with tempfile.TemporaryDirectory() as d:
            tree = _write(Path(d), self._TREE, "tree.nwk")
            m = _write(Path(d), self._PD_MATRIX, "otu_matrix.csv")
            result = faith_pd(tree, m)
            pd = dict(zip(result["samples"], result["pd"]))
            self.assertAlmostEqual(pd["S1"], 3.0)  # {A,B}: A+B+internal
            self.assertAlmostEqual(pd["S2"], 4.0)  # {A,C}: A+C+two internals
            self.assertAlmostEqual(pd["S3"], 2.0)  # {A}: A+internal

    def test_tree_layout_small_and_capped(self):
        with tempfile.TemporaryDirectory() as d:
            tree = _write(Path(d), self._TREE, "tree.nwk")
            layout = tree_layout(tree, max_tips=250)
            self.assertIsNotNone(layout)
            self.assertEqual(layout["n_tips"], 4)
            self.assertEqual(len(layout["tips"]), 4)
            self.assertIsNone(tree_layout(tree, max_tips=3))  # over cap


if __name__ == "__main__":
    unittest.main()
