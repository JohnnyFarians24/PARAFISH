"""Stage 5 (Phylo): alignment QC, branch support scaling, UniFrac."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


# Two tips per clade, so both trees have resolvable internal nodes.
FASTTREE_NEWICK = "((ASV1:0.1,ASV2:0.1)0.95:0.2,(ASV3:0.1,ASV4:0.1)0.42:0.2);"
RAXML_NEWICK = "((ASV1:0.1,ASV2:0.1)95:0.2,(ASV3:0.1,ASV4:0.1)42:0.2);"

ALIGNMENT = ">ASV1\nACGT--\n>ASV2\nACGA--\n>ASV3\nACGT-A\n>ASV4\nACGT--\n"

MATRIX = "#OTU ID,s1,s2\notu_a,10,0\notu_b,0,10\notu_c,5,5\notu_d,1,1\n"
MAPPING = "ASV_ID,Original_Header\nASV1,otu_a;size=10\nASV2,otu_b;size=10\nASV3,otu_c;size=10\nASV4,otu_d;size=2\n"


class AlignmentQcTests(unittest.TestCase):
    def test_summary_counts_sequences_length_and_gaps(self):
        from phylo.phylo_metrics import alignment_summary

        with tempfile.TemporaryDirectory() as d:
            summary = alignment_summary(_write(Path(d) / "a.fasta", ALIGNMENT))

        self.assertEqual(summary["n_sequences"], 4)
        self.assertEqual(summary["alignment_length"], 6)
        # 4 seqs x 6 cols = 24 cells; gaps: 2+2+1+2 = 7
        self.assertAlmostEqual(summary["pct_gaps"], 100 * 7 / 24, places=2)
        self.assertEqual(summary["n_reversed"], 0)

    def test_reverse_complemented_sequences_are_flagged(self):
        from phylo.phylo_metrics import alignment_summary

        with tempfile.TemporaryDirectory() as d:
            summary = alignment_summary(_write(Path(d) / "a.fasta", ">_R_ASV1\nACGT\n>ASV2\nACGT\n"))

        self.assertEqual(summary["n_reversed"], 1)
        self.assertEqual(summary["reversed_ids"], "_R_ASV1")

    def test_conservation_ignores_gaps_so_it_must_be_read_with_the_gap_fraction(self):
        """A column that is mostly gaps but invariant where present scores 1.0."""
        from phylo.phylo_metrics import alignment_columns

        with tempfile.TemporaryDirectory() as d:
            columns = alignment_columns(_write(Path(d) / "a.fasta", ">a\nA\n>b\n-\n>c\n-\n>d\n-\n"))

        self.assertEqual(columns[0]["gap_fraction"], 0.75)
        self.assertEqual(columns[0]["conservation"], 1.0)
        self.assertEqual(columns[0]["occupancy"], 1)

    def test_a_fully_variable_column_scores_low_conservation(self):
        from phylo.phylo_metrics import alignment_columns

        with tempfile.TemporaryDirectory() as d:
            columns = alignment_columns(_write(Path(d) / "a.fasta", ">a\nA\n>b\nC\n>c\nG\n>d\nT\n"))

        self.assertEqual(columns[0]["gap_fraction"], 0.0)
        self.assertEqual(columns[0]["conservation"], 0.25)

    def test_empty_alignment_yields_nothing(self):
        from phylo.phylo_metrics import alignment_columns, alignment_summary

        self.assertEqual(alignment_summary(Path("/nope.fasta")), {})
        self.assertEqual(alignment_columns(Path("/nope.fasta")), [])


class BranchSupportTests(unittest.TestCase):
    """FastTree writes 0-1, RAxML-NG writes 0-100; both must land on one scale."""

    def test_fasttree_local_support_is_scaled_to_a_percentage(self):
        from phylo.phylo_metrics import branch_support

        with tempfile.TemporaryDirectory() as d:
            rows = branch_support(_write(Path(d) / "t.nwk", FASTTREE_NEWICK))

        supports = sorted(row["support"] for row in rows)
        self.assertEqual(supports, [42.0, 95.0])

    def test_raxml_bootstrap_is_left_alone(self):
        from phylo.phylo_metrics import branch_support

        with tempfile.TemporaryDirectory() as d:
            rows = branch_support(_write(Path(d) / "t.nwk", RAXML_NEWICK))

        supports = sorted(row["support"] for row in rows)
        self.assertEqual(supports, [42.0, 95.0])

    def test_the_seventy_percent_threshold_splits_supported_clades(self):
        from phylo.phylo_metrics import branch_support

        with tempfile.TemporaryDirectory() as d:
            rows = {row["support"]: row["is_supported"] for row in branch_support(_write(Path(d) / "t.nwk", FASTTREE_NEWICK))}

        self.assertEqual(rows[95.0], 1)
        self.assertEqual(rows[42.0], 0)

    def test_a_tree_without_support_values_yields_nothing(self):
        from phylo.phylo_metrics import branch_support

        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(branch_support(_write(Path(d) / "t.nwk", "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);")), [])

    def test_a_missing_tree_yields_nothing(self):
        from phylo.phylo_metrics import branch_support

        self.assertEqual(branch_support(Path("/nope.nwk")), [])


class UnifracTests(unittest.TestCase):
    def _fixture(self, root: Path):
        return (
            _write(root / "t.nwk", FASTTREE_NEWICK),
            _write(root / "m.csv", MATRIX),
            _write(root / "map.csv", MAPPING),
        )

    def test_weighted_and_unweighted_differ_on_the_same_data(self):
        """Unweighted counts rare lineages as much as dominant ones; weighted does not."""
        from phylo.phylo_metrics import unifrac_distances

        with tempfile.TemporaryDirectory() as d:
            tree, matrix, mapping = self._fixture(Path(d))
            weighted = unifrac_distances(tree, matrix, mapping, weighted=True)
            unweighted = unifrac_distances(tree, matrix, mapping, weighted=False)

        self.assertEqual(weighted["samples"], ["s1", "s2"])
        self.assertEqual(weighted["matrix"][0][0], 0.0)
        self.assertGreater(weighted["matrix"][0][1], 0.0)
        self.assertNotAlmostEqual(weighted["matrix"][0][1], unweighted["matrix"][0][1], places=3)

    def test_the_distance_matrix_is_symmetric_with_a_zero_diagonal(self):
        from phylo.phylo_metrics import unifrac_distances

        with tempfile.TemporaryDirectory() as d:
            tree, matrix, mapping = self._fixture(Path(d))
            result = unifrac_distances(tree, matrix, mapping)

        values = result["matrix"]
        self.assertEqual(values[0][0], 0.0)
        self.assertEqual(values[1][1], 0.0)
        self.assertAlmostEqual(values[0][1], values[1][0], places=6)

    def test_without_the_id_mapping_no_tip_matches_a_feature(self):
        """Tips are ASV1..N but features are the original OTU ids."""
        from phylo.phylo_metrics import unifrac_distances

        with tempfile.TemporaryDirectory() as d:
            tree, matrix, _ = self._fixture(Path(d))
            self.assertEqual(unifrac_distances(tree, matrix, None), {})

    def test_a_single_sample_has_no_between_sample_distance(self):
        from phylo.phylo_metrics import unifrac_distances

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            tree = _write(root / "t.nwk", FASTTREE_NEWICK)
            matrix = _write(root / "m.csv", "#OTU ID,s1\notu_a,10\notu_b,5\n")
            mapping = _write(root / "map.csv", MAPPING)
            self.assertEqual(unifrac_distances(tree, matrix, mapping), {})

    def test_features_mapped_to_tips_absent_from_the_tree_are_dropped(self):
        """A stale mapping must not silently become a matrix of zeros."""
        from phylo.phylo_metrics import unifrac_distances

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            tree = _write(root / "t.nwk", FASTTREE_NEWICK)
            matrix = _write(root / "m.csv", MATRIX)
            # otu_d points at ASV99, which this tree does not contain.
            stale = _write(root / "map.csv", MAPPING.replace("ASV4,otu_d;size=2", "ASV99,otu_d;size=2"))
            result = unifrac_distances(tree, matrix, stale)

        self.assertEqual(result["samples"], ["s1", "s2"])
        self.assertGreater(result["matrix"][0][1], 0.0)

    def test_a_mapping_that_matches_no_tip_yields_nothing(self):
        from phylo.phylo_metrics import unifrac_distances

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            tree = _write(root / "t.nwk", FASTTREE_NEWICK)
            matrix = _write(root / "m.csv", MATRIX)
            bogus = _write(root / "map.csv", "ASV_ID,Original_Header\nZZZ1,otu_a\nZZZ2,otu_b\nZZZ3,otu_c\nZZZ4,otu_d\n")
            self.assertEqual(unifrac_distances(tree, matrix, bogus), {})


class PhyloPersistenceTests(unittest.TestCase):
    def test_build_writes_every_table(self):
        from phylo.phylo_metrics import build_phylo_metrics

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            files = build_phylo_metrics(
                root,
                _write(root / "t.nwk", FASTTREE_NEWICK),
                _write(root / "a.fasta", ALIGNMENT),
                _write(root / "m.csv", MATRIX),
                _write(root / "map.csv", MAPPING),
            )

            self.assertIn("alignment_summary_csv", files)
            self.assertIn("alignment_columns_csv", files)
            self.assertIn("branch_support_csv", files)
            self.assertIn("weighted_unifrac_csv", files)
            self.assertIn("unweighted_unifrac_csv", files)

            with open(files["weighted_unifrac_csv"]) as handle:
                header = next(csv.reader(handle))
            self.assertEqual(header, ["sample", "s1", "s2"])

    def test_alignment_only_run_skips_the_tree_tables(self):
        from phylo.phylo_metrics import build_phylo_metrics

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            files = build_phylo_metrics(root, aligned_fasta=_write(root / "a.fasta", ALIGNMENT))

        self.assertIn("alignment_summary_csv", files)
        self.assertNotIn("branch_support_csv", files)


class SupportTreeSelectionTests(unittest.TestCase):
    """RAxML's bestTree is bare; branch support must read the *.support file."""

    def test_the_on_disk_support_file_beats_a_bare_best_tree(self):
        from phylo.phylo_metrics import support_tree_path

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            best = _write(root / "raxml" / "raxml.raxml.bestTree", RAXML_NEWICK.replace("95", "").replace("42", ""))
            support = _write(root / "raxml" / "raxml.raxml.support", RAXML_NEWICK)
            # The record points support_tree at a name RAxML did not actually write.
            files = {"support_tree": str(root / "raxml" / "raxml.raxml.supportFBP"), "tree_nwk": str(best), "best_tree": str(best)}
            chosen = support_tree_path(files, root)

        self.assertEqual(Path(chosen), support)

    def test_a_valid_recorded_support_tree_is_used_directly(self):
        from phylo.phylo_metrics import support_tree_path

        with tempfile.TemporaryDirectory() as d:
            support = _write(Path(d) / "s.nwk", RAXML_NEWICK)
            self.assertEqual(support_tree_path({"support_tree": str(support)}), support)

    def test_fasttree_tree_is_used_when_there_is_no_raxml_support(self):
        from phylo.phylo_metrics import support_tree_path

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            tree = _write(root / "fasttree" / "tree.nwk", FASTTREE_NEWICK)
            self.assertEqual(support_tree_path({"tree_nwk": str(tree)}, root), tree)


class PhyloRunnerFileTests(unittest.TestCase):
    def test_the_alignment_artifacts_reach_the_final_step_output(self):
        """run_phylo used to return only the tree builder's files."""
        from core.models import StepOutput
        from phylo.phylo_runner import _with_alignment_files

        mafft = StepOutput(step_name="mafft", success=True, files={"aligned_fasta": Path("/x/aligned.fasta"), "reversed_asvs_csv": Path("/x/map.csv")})
        tree = StepOutput(step_name="fasttree", success=True, files={"tree_nwk": Path("/x/tree.nwk")})
        merged = _with_alignment_files(tree, mafft)

        self.assertEqual(merged.files["aligned_fasta"], Path("/x/aligned.fasta"))
        self.assertEqual(merged.files["reversed_asvs_csv"], Path("/x/map.csv"))
        self.assertEqual(merged.files["tree_nwk"], Path("/x/tree.nwk"))

    def test_existing_keys_are_not_overwritten(self):
        from core.models import StepOutput
        from phylo.phylo_runner import _with_alignment_files

        mafft = StepOutput(step_name="mafft", success=True, files={"aligned_fasta": Path("/x/a.fasta")})
        tree = StepOutput(step_name="raxml", success=True, files={"aligned_fasta": Path("/y/own.fasta")})

        self.assertEqual(_with_alignment_files(tree, mafft).files["aligned_fasta"], Path("/y/own.fasta"))


class PhyloFigureTests(unittest.TestCase):
    def _layout(self, tree_path: Path):
        from phylo.phylo_diversity import tree_layout

        return tree_layout(tree_path)

    def test_tree_layout_exposes_internal_node_support(self):
        with tempfile.TemporaryDirectory() as d:
            layout = self._layout(_write(Path(d) / "t.nwk", FASTTREE_NEWICK))

        supports = sorted(round(node[2], 2) for node in layout["nodes"])
        self.assertEqual(supports, [0.42, 0.95])

    def test_annotated_tree_splits_nodes_by_the_support_threshold(self):
        from ui.plots import annotated_tree_figure

        with tempfile.TemporaryDirectory() as d:
            layout = self._layout(_write(Path(d) / "t.nwk", FASTTREE_NEWICK))
        names = [trace.name for trace in annotated_tree_figure(layout).data if trace.name]
        self.assertIn("support >= 70%", names)
        self.assertIn("support < 70%", names)

    def test_abundance_bars_only_appear_when_abundance_is_known(self):
        from ui.plots import annotated_tree_figure

        with tempfile.TemporaryDirectory() as d:
            layout = self._layout(_write(Path(d) / "t.nwk", FASTTREE_NEWICK))

        without = annotated_tree_figure(layout)
        with_bars = annotated_tree_figure(layout, abundance={"ASV1": 10, "ASV2": 5, "ASV3": 1, "ASV4": 0})
        self.assertNotIn("bar", [trace.type for trace in without.data])
        self.assertIn("bar", [trace.type for trace in with_bars.data])

    def test_empty_layout_renders_nothing(self):
        from ui.plots import annotated_tree_figure

        self.assertIsNone(annotated_tree_figure({}))
        self.assertIsNone(annotated_tree_figure(None))


if __name__ == "__main__":
    unittest.main()
