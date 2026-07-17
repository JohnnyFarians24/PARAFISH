"""Stage 3 (Analysis): chimeras, rarefaction, alpha diversity, UpSet data."""

from __future__ import annotations

import csv
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _matrix_csv(path: Path, samples: list[str], rows: dict[str, list[int]]) -> Path:
    lines = ["#OTU ID," + ",".join(samples)]
    lines += [f"{feature}," + ",".join(str(v) for v in counts) for feature, counts in rows.items()]
    return _write(path, "\n".join(lines) + "\n")


class ChimeraStatsTests(unittest.TestCase):
    """The variant share and the read share answer different questions."""

    def _fixture(self, root: Path):
        # Two centroids: c1 carries a huge cluster, c2 (the chimera) a tiny one.
        _write(root / "otus_raw.fasta", ">c1;size=10\nACGT\n>c2;size=5\nTTTT\n")
        _write(root / "chimeras.fasta", ">c2;size=5\nTTTT\n")
        _write(
            root / "clusters.uc",
            "S\t0\t4\t*\t*\t*\t*\t*\tc1;size=10\t*\n"
            "H\t0\t4\t99.0\t+\t0\t0\t4M\tm1;size=990\tc1;size=10\n"
            "S\t1\t4\t*\t*\t*\t*\t*\tc2;size=5\t*\n"
            "H\t1\t4\t99.0\t+\t0\t0\t4M\tm2;size=5\tc2;size=5\n",
        )
        return root

    def test_variant_share_is_independent_of_read_share(self):
        from analysis.analysis_metrics import chimera_stats

        with tempfile.TemporaryDirectory() as d:
            root = self._fixture(Path(d))
            stats = chimera_stats(root / "chimeras.fasta", root / "otus_raw.fasta", root / "clusters.uc")

        # Half the variants are chimeric...
        self.assertEqual(stats["chimeric_variants"], 1)
        self.assertEqual(stats["total_variants"], 2)
        self.assertEqual(stats["pct_chimeric_variants"], 50.0)
        # ...but they carry only 10 of 1010 clustered reads.
        self.assertEqual(stats["chimeric_reads"], 10)
        self.assertEqual(stats["total_reads"], 1010)
        self.assertAlmostEqual(stats["pct_chimeric_reads"], 0.99, places=2)
        self.assertEqual(stats["read_basis"], "cluster")

    def test_total_reads_is_the_canonical_mapping_total_when_a_matrix_exists(self):
        """The structural bug: chimera_stats counted a different population.

        sum(cluster_read_counts) excludes every read --minuniquesize dropped at
        dereplication, so this file reported 147 010 while otu_table, taxonomic_
        resolution and the taxa_tables all reported 157 527. The abundance table is
        canonical; chimera_stats must quote it.
        """
        from analysis.analysis_metrics import chimera_stats

        with tempfile.TemporaryDirectory() as d:
            root = self._fixture(Path(d))
            matrix = root / "otu_matrix.csv"
            # Deliberately != 1010: the mapping sees reads the clustering never did.
            matrix.write_text("#OTU ID,s1,s2\nOTU1,900,300\n", encoding="utf-8")
            stats = chimera_stats(root / "chimeras.fasta", root / "otus_raw.fasta", root / "clusters.uc", matrix)

        self.assertEqual(stats["total_reads"], 1200)
        self.assertEqual(stats["read_basis"], "otu_table")
        # Chimeric reads are absent from otu_table by construction (chimeras are
        # removed before mapping), so the denominator adds them back: 10/(1200+10).
        self.assertEqual(stats["chimeric_reads"], 10)
        self.assertAlmostEqual(stats["pct_chimeric_reads"], 0.826, places=3)

    def test_cluster_fallback_does_not_double_count_the_chimeric_reads(self):
        """sum(cluster_read_counts) already contains them; the otu_table total does not.
        The denominator therefore depends on the basis, and conflating the two inflates
        or deflates the rate."""
        from analysis.analysis_metrics import chimera_stats

        with tempfile.TemporaryDirectory() as d:
            root = self._fixture(Path(d))
            stats = chimera_stats(root / "chimeras.fasta", root / "otus_raw.fasta", root / "clusters.uc")

        self.assertEqual(stats["read_basis"], "cluster")
        self.assertAlmostEqual(stats["pct_chimeric_reads"], 0.99, places=2)  # 10/1010, not 10/1020

    def test_without_the_cluster_map_the_denominator_is_the_centroids_only(self):
        """The centroid's own size= is its dereplicated abundance, not its cluster's."""
        from analysis.analysis_metrics import chimera_stats

        with tempfile.TemporaryDirectory() as d:
            root = self._fixture(Path(d))
            stats = chimera_stats(root / "chimeras.fasta", root / "otus_raw.fasta")

        self.assertEqual(stats["total_reads"], 15)  # 10 + 5, not 1010
        self.assertEqual(stats["read_basis"], "centroid")

    def test_cluster_read_counts_sum_to_every_dereplicated_read(self):
        from analysis.analysis_metrics import cluster_read_counts

        with tempfile.TemporaryDirectory() as d:
            root = self._fixture(Path(d))
            counts = cluster_read_counts(root / "clusters.uc")

        self.assertEqual(counts, {"c1;size=10": 1000, "c2;size=5": 10})

    def test_no_chimeras_reports_zero(self):
        from analysis.analysis_metrics import chimera_stats

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _write(root / "otus_raw.fasta", ">c1;size=10\nACGT\n")
            _write(root / "chimeras.fasta", "")
            stats = chimera_stats(root / "chimeras.fasta", root / "otus_raw.fasta")

        self.assertEqual(stats["pct_chimeric_variants"], 0.0)
        self.assertEqual(stats["pct_chimeric_reads"], 0.0)


class AlphaDiversityTests(unittest.TestCase):
    def test_hill_numbers_agree_with_the_classic_indices(self):
        """hill_q0 == richness and hill_q1 == exp(Shannon) only if both use ln."""
        from analysis.analysis_metrics import alpha_diversity_row

        counts = np.array([5, 0, 3, 10, 1, 1, 2])
        row = alpha_diversity_row("s1", counts)

        self.assertEqual(row["observed"], 6)
        self.assertEqual(row["hill_q0"], 6.0)
        self.assertAlmostEqual(row["hill_q1"], math.exp(row["shannon"]), places=2)

    def test_chao1_uses_the_bias_corrected_estimator(self):
        from analysis.analysis_metrics import alpha_diversity_row

        # 6 observed, 2 singletons, 1 doubleton -> 6 + (2*1)/(2*(1+1)) = 6.5
        row = alpha_diversity_row("s1", np.array([5, 0, 3, 10, 1, 1, 2]))
        self.assertEqual(row["chao1"], 6.5)

    def test_a_single_feature_has_zero_evenness_and_zero_shannon(self):
        from analysis.analysis_metrics import alpha_diversity_row

        row = alpha_diversity_row("s1", np.array([100, 0, 0]))
        self.assertEqual(row["observed"], 1)
        self.assertEqual(row["shannon"], 0.0)
        self.assertEqual(row["pielou"], 0.0)

    def test_empty_sample_does_not_divide_by_zero(self):
        from analysis.analysis_metrics import alpha_diversity_row

        row = alpha_diversity_row("s1", np.array([0, 0]))
        self.assertEqual(row["reads"], 0)
        self.assertEqual(row["shannon"], 0.0)


class RarefactionTests(unittest.TestCase):
    def test_rarefy_draws_exactly_the_requested_depth(self):
        from analysis.analysis_metrics import rarefy

        rng = np.random.default_rng(1)
        counts = np.array([50, 30, 20])
        drawn = rarefy(counts, 40, rng)
        self.assertEqual(int(drawn.sum()), 40)
        self.assertTrue((drawn <= counts).all())

    def test_rarefy_is_a_noop_when_the_sample_is_shallower_than_the_depth(self):
        from analysis.analysis_metrics import rarefy

        counts = np.array([5, 3])
        drawn = rarefy(counts, 100, np.random.default_rng(1))
        self.assertTrue((drawn == counts).all())

    def test_every_sample_is_rarefied_to_the_shallowest_depth(self):
        from analysis.analysis_metrics import build_analysis_metrics

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            matrix = _matrix_csv(root / "otu_matrix.csv", ["deep", "shallow"], {"f1": [800, 40], "f2": [200, 60]})
            files = build_analysis_metrics(root, matrix)
            depth = list(csv.DictReader(open(files["rarefaction_depth_csv"])))[0]
            alpha = list(csv.DictReader(open(files["alpha_diversity_csv"])))

        self.assertEqual(depth["depth"], "100")
        self.assertEqual({row["reads"] for row in alpha}, {"100"})

    def test_empty_samples_are_dropped_not_used_as_the_depth(self):
        from analysis.analysis_metrics import build_analysis_metrics

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            matrix = _matrix_csv(root / "otu_matrix.csv", ["ok", "empty"], {"f1": [50, 0], "f2": [50, 0]})
            files = build_analysis_metrics(root, matrix)
            depth = list(csv.DictReader(open(files["rarefaction_depth_csv"])))[0]
            alpha = list(csv.DictReader(open(files["alpha_diversity_csv"])))

        self.assertEqual(depth["depth"], "100")
        self.assertEqual(depth["n_samples_dropped"], "1")
        self.assertEqual([row["sample"] for row in alpha], ["ok"])

    def test_rarefaction_is_reproducible_for_a_given_seed(self):
        from analysis.analysis_metrics import build_analysis_metrics

        def _run(root: Path):
            matrix = _matrix_csv(root / "m.csv", ["a", "b"], {"f1": [80, 40], "f2": [20, 60]})
            files = build_analysis_metrics(root, matrix, seed=7)
            return Path(files["alpha_diversity_csv"]).read_text()

        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            self.assertEqual(_run(Path(d1)), _run(Path(d2)))


class MatrixAndSetTests(unittest.TestCase):
    def test_matrix_loads_as_samples_by_features(self):
        from analysis.analysis_metrics import load_matrix

        with tempfile.TemporaryDirectory() as d:
            matrix = _matrix_csv(Path(d) / "m.csv", ["s1", "s2"], {"f1": [1, 2], "f2": [3, 4]})
            features, samples, counts = load_matrix(matrix)

        self.assertEqual(features, ["f1", "f2"])
        self.assertEqual(samples, ["s1", "s2"])
        self.assertEqual(counts.shape, (2, 2))
        self.assertEqual(counts[0].tolist(), [1.0, 3.0])  # sample s1 across features

    def test_intersections_partition_the_features(self):
        from analysis.analysis_metrics import feature_intersections, load_matrix

        with tempfile.TemporaryDirectory() as d:
            matrix = _matrix_csv(
                Path(d) / "m.csv", ["a", "b"], {"only_a": [5, 0], "only_b": [0, 5], "shared": [1, 1], "shared2": [2, 3]}
            )
            _, samples, counts = load_matrix(matrix)
            rows = feature_intersections(counts, samples)

        by_combination = {row["combination"]: row["n_features"] for row in rows}
        self.assertEqual(by_combination, {"a": 1, "b": 1, "a|b": 2})
        self.assertEqual(sum(by_combination.values()), 4)

    def test_core_features_use_the_prevalence_threshold(self):
        from analysis.analysis_metrics import feature_prevalence, load_matrix

        with tempfile.TemporaryDirectory() as d:
            matrix = _matrix_csv(Path(d) / "m.csv", ["a", "b"], {"everywhere": [1, 1], "rare": [1, 0]})
            features, _, counts = load_matrix(matrix)
            rows = {row["feature"]: row for row in feature_prevalence(counts, features, core_fraction=0.8)}

        self.assertEqual(rows["everywhere"]["is_core"], 1)
        self.assertEqual(rows["rare"]["is_core"], 0)
        self.assertEqual(rows["rare"]["prevalence"], 0.5)

    def test_rank_abundance_is_relative_and_descending(self):
        from analysis.analysis_metrics import load_matrix, rank_abundance

        with tempfile.TemporaryDirectory() as d:
            matrix = _matrix_csv(Path(d) / "m.csv", ["a"], {"f1": [10], "f2": [30], "f3": [0]})
            _, samples, counts = load_matrix(matrix)
            rows = rank_abundance(counts, samples)

        self.assertEqual([row["rank"] for row in rows], [1, 2])
        self.assertEqual(rows[0]["relative_abundance"], 0.75)
        self.assertGreater(rows[0]["relative_abundance"], rows[1]["relative_abundance"])


class FastaHelperTests(unittest.TestCase):
    def test_sizes_default_to_one_when_absent(self):
        from analysis.analysis_metrics import fasta_size_total

        with tempfile.TemporaryDirectory() as d:
            path = _write(Path(d) / "x.fasta", ">a;size=7\nAC\n>b\nGT\n")
            self.assertEqual(fasta_size_total(path), (2, 8))

    def test_lengths_span_wrapped_sequences(self):
        from analysis.analysis_metrics import fasta_lengths

        with tempfile.TemporaryDirectory() as d:
            path = _write(Path(d) / "x.fasta", ">a\nACGT\nACGT\n>b\nAC\n")
            self.assertEqual(dict(fasta_lengths(path)), {8: 1, 2: 1})

    def test_missing_fasta_is_empty(self):
        from analysis.analysis_metrics import fasta_lengths, fasta_size_total

        self.assertEqual(fasta_size_total(Path("/nope.fasta")), (0, 0))
        self.assertEqual(dict(fasta_lengths(Path("/nope.fasta"))), {})


if __name__ == "__main__":
    unittest.main()
