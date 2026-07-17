"""Ordination: compositional correctness, and refusing to draw a misleading plot."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from analysis.ordination import (
    LOW_CONFIDENCE_SAMPLES,
    MIN_SAMPLES,
    OrdinationError,
    clr,
    filter_features,
    ordinate,
    pca,
    permanova,
    top_loadings,
    write_coordinates,
    write_distance_matrix,
)


def _two_groups(seed: int = 0, per_group: int = 3):
    """Two communities with opposite dominance — a separation that must be recovered."""
    rng = np.random.default_rng(seed)
    a = rng.multinomial(10_000, [0.40, 0.30, 0.15, 0.10, 0.05], size=per_group)
    b = rng.multinomial(10_000, [0.05, 0.10, 0.15, 0.30, 0.40], size=per_group)
    counts = np.vstack([a, b]).astype(float)
    samples = [f"A{i}" for i in range(per_group)] + [f"B{i}" for i in range(per_group)]
    features = [f"OTU{i}" for i in range(1, 6)]
    return counts, samples, features


class GuardTests(unittest.TestCase):
    def test_two_samples_are_refused_not_plotted(self):
        """With 2 samples there is one axis and the 'separation' is an artefact of n."""
        with self.assertRaises(OrdinationError) as caught:
            ordinate(np.array([[10.0, 5.0, 1.0], [3.0, 8.0, 2.0]]), ["A", "B"], ["F1", "F2", "F3"])
        self.assertIn(str(MIN_SAMPLES), str(caught.exception))

    def test_three_samples_plot_but_are_flagged_low_confidence(self):
        counts, samples, features = _two_groups(per_group=2)
        counts, samples = counts[:3], samples[:3]
        result = ordinate(counts, samples, features)
        self.assertTrue(result["low_confidence"])
        self.assertLess(result["n_samples"], LOW_CONFIDENCE_SAMPLES)

    def test_too_few_features_is_refused(self):
        counts = np.array([[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]])
        # Every feature is in exactly one sample, so the prevalence filter empties it.
        with self.assertRaises(OrdinationError):
            ordinate(counts, ["A", "B", "C"], ["F1", "F2", "F3"], min_prevalence=2)


class CompositionalCorrectnessTests(unittest.TestCase):
    def test_clr_is_invariant_to_sequencing_depth(self):
        """The whole reason plain PCA is wrong here.

        Two samples with the SAME community but different sequencing depth must land in
        the same place. On raw counts they do not — PC1 becomes 'how deep did we
        sequence', which is a library artefact, not biology.
        """
        composition = np.array([[40.0, 30, 15, 10, 5]])
        shallow = composition * 100
        deep = composition * 5_000

        transformed = clr(np.vstack([shallow, deep]))
        np.testing.assert_allclose(transformed[0], transformed[1], rtol=1e-2)

    def test_raw_pca_captures_depth_but_clr_pca_does_not(self):
        counts, samples, features = _two_groups()
        deep = counts.copy()
        deep[3:] *= 50  # same composition, 50x the reads

        # On raw counts, the dominant axis is depth: the two blocks separate on PC1
        # purely because one was sequenced deeper.
        raw_coords, _ratio, _v = pca(deep)
        depth_gap = abs(raw_coords[3:, 0].mean() - raw_coords[:3, 0].mean())
        spread_within = float(np.std(raw_coords[:3, 0]))
        self.assertGreater(depth_gap, spread_within * 10)

        # CLR is scale-invariant, so a 50x depth change barely moves the samples. It is
        # *approximately*, not exactly, invariant — see test_pseudocount_breaks_exact_
        # scale_invariance below — so the tolerance is loose but still far tighter than
        # the depth artefact it replaces.
        plain = ordinate(counts, samples, features, metric="aitchison")
        scaled = ordinate(deep, samples, features, metric="aitchison")
        np.testing.assert_allclose(
            np.abs(plain["coordinates"]), np.abs(scaled["coordinates"]), rtol=0.01
        )

    def test_pseudocount_breaks_exact_scale_invariance(self):
        """Documents a real property, so nobody later "fixes" it by tightening a test.

        ``clr(x)`` is exactly scale-invariant. ``clr(x + 0.5)`` is not: adding a constant
        0.5 means something different at a depth of 100 than at 5,000. The distortion is
        tiny (well under 1% here, for a 50x depth swing) and orders of magnitude below the
        biological signal — but it is nonzero, and it is the price of a fixed pseudocount.
        Bayesian multiplicative replacement (zCompositions) is what removes it.
        """
        composition = np.array([[40.0, 30, 15, 10, 5]])
        shallow = clr(composition * 100)
        deep = clr(composition * 5_000)

        difference = float(np.abs(shallow - deep).max())
        self.assertGreater(difference, 0.0)     # not exact...
        self.assertLess(difference, 0.05)       # ...but negligible in practice

    def test_known_groups_separate_on_pc1(self):
        counts, samples, features = _two_groups()
        for metric in ("aitchison", "braycurtis"):
            result = ordinate(counts, samples, features, metric=metric)
            pc1 = result["coordinates"][:, 0]
            # The two blocks must sit on opposite sides of the origin.
            self.assertLess(pc1[:3].mean() * pc1[3:].mean(), 0, metric)

    def test_explained_variance_is_reported_and_sums_sanely(self):
        counts, samples, features = _two_groups()
        result = ordinate(counts, samples, features)
        self.assertGreater(result["explained"][0], 0)
        self.assertLessEqual(result["explained"][0], 1.0 + 1e-9)

    def test_degenerate_distance_matrix_is_refused(self):
        """Jaccard sees only presence: if every feature is in every sample, all the
        distances are zero and skbio hands back NaN axes. Plotting those is nonsense."""
        counts = np.array([[5.0, 3, 2], [7.0, 1, 9], [2.0, 8, 4]])
        with self.assertRaises(OrdinationError):
            ordinate(counts, ["A", "B", "C"], ["F1", "F2", "F3"], metric="jaccard", min_prevalence=1)


class PrevalenceFilterTests(unittest.TestCase):
    def test_singletons_are_dropped(self):
        counts = np.array([[5.0, 1.0, 0.0], [3.0, 0.0, 0.0], [2.0, 4.0, 0.0]])
        filtered, features = filter_features(counts, ["F1", "F2", "F3"], min_prevalence=2)
        self.assertEqual(features, ["F1", "F2"])  # F3 absent, F2 in 2 samples
        self.assertEqual(filtered.shape[1], 2)


class LoadingsAndStatsTests(unittest.TestCase):
    def test_loadings_are_named_by_taxon_when_known(self):
        counts, samples, features = _two_groups()
        result = ordinate(counts, samples, features)
        loadings = top_loadings(result, {"OTU1": "Diplostomum spathaceum"}, top_n=5)
        labels = {entry["feature"]: entry["label"] for entry in loadings}
        self.assertEqual(labels["OTU1"], "Diplostomum spathaceum")
        self.assertEqual(labels["OTU2"], "OTU2")  # no taxon: falls back to the id

    def test_permanova_needs_three_replicates_per_group(self):
        counts, samples, features = _two_groups(per_group=2)
        result = ordinate(counts, samples, features)
        groups = {sample: sample[0] for sample in samples}  # 2 + 2
        self.assertIsNone(permanova(result, groups))

    def test_permanova_reports_its_p_value_floor(self):
        """With 3 v 3 the smallest attainable p is 1/10 — saying 'p = 0.12' without that
        context invites the reader to conclude there is no effect."""
        counts, samples, features = _two_groups(per_group=3)
        result = ordinate(counts, samples, features)
        stats = permanova(result, {sample: sample[0] for sample in samples})
        self.assertIsNotNone(stats)
        self.assertAlmostEqual(stats["p_floor"], 0.1, places=6)
        self.assertGreater(stats["pseudo_f"], 0)


class ExportTests(unittest.TestCase):
    def test_exported_coordinates_match_the_plotted_ones(self):
        counts, samples, features = _two_groups()
        result = ordinate(counts, samples, features)
        with tempfile.TemporaryDirectory() as temp_dir:
            path = write_coordinates(result, Path(temp_dir) / "coords.tsv")
            rows = [line.split("\t") for line in path.read_text().splitlines() if line and not line.startswith("#")]

        body = [row for row in rows[1:] if len(row) == 3]
        self.assertEqual(len(body), len(samples))
        for index, row in enumerate(body):
            self.assertEqual(row[0], samples[index])
            self.assertAlmostEqual(float(row[1]), result["coordinates"][index, 0], places=5)
            self.assertAlmostEqual(float(row[2]), result["coordinates"][index, 1], places=5)

    def test_distance_matrix_export_is_square_and_symmetric(self):
        counts, samples, features = _two_groups()
        result = ordinate(counts, samples, features, metric="braycurtis")
        with tempfile.TemporaryDirectory() as temp_dir:
            path = write_distance_matrix(result, Path(temp_dir) / "dist.tsv")
            rows = [line.split("\t") for line in path.read_text().splitlines() if line]

        self.assertEqual(rows[0][1:], samples)
        values = np.array([[float(cell) for cell in row[1:]] for row in rows[1:]])
        self.assertEqual(values.shape, (len(samples), len(samples)))
        np.testing.assert_allclose(values, values.T, atol=1e-9)
        np.testing.assert_allclose(np.diag(values), 0, atol=1e-9)


if __name__ == "__main__":
    unittest.main()
