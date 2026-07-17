"""Stage 4 (Databases): resolution, confidence, decontamination, prevalence."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


CONSENSUS_HEADER = (
    "ASV_ID,Consensus_Kingdom,Consensus_Phylum,Consensus_Class,Consensus_Order,"
    "Consensus_Family,Consensus_Genus,Consensus_Species\n"
)


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _consensus(path: Path) -> Path:
    return _write(
        path,
        CONSENSUS_HEADER
        + "to_species,Eukaryota,P,C,O,F,Diplostomum,Diplostomum spathaceum\n"
        + "to_genus,Bacteria,P,C,O,F,Ralstonia,\n"
        + "to_kingdom,Eukaryota,,,,,,\n",
    )


def _matrix(path: Path) -> Path:
    return _write(
        path,
        "#OTU ID,ctrl1,ctrl2,s1,s2,s3,s4\n"
        "to_species,0,0,100,120,90,80\n"
        "to_genus,10,8,5,0,0,0\n"
        "to_kingdom,0,0,3,0,0,0\n",
    )


class ReliabilityTests(unittest.TestCase):
    """Identity is graded against the rank it is being used to support."""

    def test_high_identity_and_coverage_is_high_confidence(self):
        from databases.taxonomy_runner import assignment_reliability

        self.assertEqual(assignment_reliability({"similarity": "99", "coverage": "100"}, "species"), "high")

    def test_high_identity_over_a_short_alignment_is_only_moderate(self):
        from databases.taxonomy_runner import assignment_reliability

        self.assertEqual(assignment_reliability({"similarity": "99", "coverage": "40"}, "species"), "moderate")

    def test_the_same_identity_supports_a_family_but_not_a_species(self):
        from databases.taxonomy_runner import assignment_reliability

        record = {"similarity": "91", "coverage": "100"}
        self.assertEqual(assignment_reliability(record, "species"), "low")
        self.assertEqual(assignment_reliability(record, "family"), "high")

    def test_sintax_is_graded_on_its_bootstrap(self):
        from databases.taxonomy_runner import assignment_reliability

        self.assertEqual(assignment_reliability({"confidence": "0.95"}, "genus"), "high")
        self.assertEqual(assignment_reliability({"confidence": "0.54"}, "species"), "low")

    def test_a_source_reporting_no_evidence_cannot_be_graded(self):
        from databases.taxonomy_runner import assignment_reliability

        self.assertEqual(assignment_reliability({}, "genus"), "unknown")

    def test_unassigned_records_are_not_graded(self):
        from databases.taxonomy_runner import assignment_reliability

        self.assertEqual(assignment_reliability({"similarity": "99"}, "unassigned"), "unassigned")

    def test_sintax_records_carry_no_identity_or_coverage(self):
        """A k-mer classifier does not align, so inventing coverage would be a lie."""
        from databases.taxonomy_runner import _collect_sintax_records

        with tempfile.TemporaryDirectory() as d:
            path = _write(Path(d) / "s.csv", "asv1\td:Eukaryota(1.00),g:Salmo(0.91)\t+\td:Eukaryota,g:Salmo\n")
            record = _collect_sintax_records("PR2", path)["asv1"][0]

        self.assertEqual(record["similarity"], "")
        self.assertEqual(record["coverage"], "")
        self.assertEqual(record["confidence"], "0.91")


class ResolutionTests(unittest.TestCase):
    def test_features_and_reads_are_counted_separately(self):
        from databases.database_metrics import taxonomic_resolution

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            rows = {row["rank"]: row for row in taxonomic_resolution(_consensus(root / "c.csv"), _matrix(root / "m.csv"))}

        # One feature each at species/genus/kingdom, but the reads are wildly skewed.
        self.assertEqual(rows["species"]["n_features"], 1)
        self.assertEqual(rows["genus"]["n_features"], 1)
        self.assertEqual(rows["kingdom"]["n_features"], 1)
        self.assertEqual(rows["species"]["n_reads"], 390)
        self.assertEqual(rows["genus"]["n_reads"], 23)
        self.assertAlmostEqual(rows["species"]["pct_features"], 33.333, places=2)
        self.assertGreater(rows["species"]["pct_reads"], 90.0)

    def test_deepest_rank_stops_at_the_last_filled_column(self):
        from databases.database_metrics import deepest_consensus_rank

        self.assertEqual(deepest_consensus_rank({"Consensus_Kingdom": "Eukaryota", "Consensus_Genus": "Salmo"}), "genus")
        self.assertEqual(deepest_consensus_rank({"Consensus_Kingdom": "Eukaryota"}), "kingdom")
        self.assertEqual(deepest_consensus_rank({"Consensus_Kingdom": "unclassified"}), "unassigned")

    def test_resolution_without_a_matrix_still_counts_features(self):
        from databases.database_metrics import taxonomic_resolution

        with tempfile.TemporaryDirectory() as d:
            rows = {row["rank"]: row for row in taxonomic_resolution(_consensus(Path(d) / "c.csv"))}

        self.assertEqual(rows["species"]["n_features"], 1)
        self.assertEqual(rows["species"]["n_reads"], 0)


class WilsonIntervalTests(unittest.TestCase):
    def test_matches_the_published_interval(self):
        from databases.database_metrics import wilson_interval

        low, high = wilson_interval(3, 5)
        self.assertAlmostEqual(low, 0.2306, places=3)
        self.assertAlmostEqual(high, 0.8824, places=3)

    def test_never_leaves_the_unit_interval_at_the_extremes(self):
        """The normal approximation returns negative lower bounds here."""
        from databases.database_metrics import wilson_interval

        low, high = wilson_interval(0, 5)
        self.assertEqual(low, 0.0)
        self.assertLess(high, 1.0)
        low, high = wilson_interval(5, 5)
        self.assertGreater(low, 0.0)
        self.assertEqual(high, 1.0)

    def test_zero_trials_does_not_divide_by_zero(self):
        from databases.database_metrics import wilson_interval

        self.assertEqual(wilson_interval(0, 0), (0.0, 0.0))


class DecontamTests(unittest.TestCase):
    def test_no_controls_means_no_calls_rather_than_guessing(self):
        from databases.database_metrics import decontam_prevalence

        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(decontam_prevalence(_matrix(Path(d) / "m.csv"), []), [])

    def test_a_feature_confined_to_the_controls_scores_lowest(self):
        from databases.database_metrics import decontam_prevalence

        with tempfile.TemporaryDirectory() as d:
            matrix = _write(
                Path(d) / "m.csv",
                "#OTU ID,c1,c2,c3,s1,s2,s3\n"
                "reagent,5,5,5,0,0,0\n"
                "biological,0,0,0,50,60,70\n",
            )
            rows = {row["feature"]: row for row in decontam_prevalence(matrix, ["c1", "c2", "c3"])}

        self.assertEqual(rows["reagent"]["prevalence_controls"], 1.0)
        self.assertEqual(rows["reagent"]["prevalence_samples"], 0.0)
        self.assertLess(rows["reagent"]["score"], rows["biological"]["score"])
        self.assertEqual(rows["biological"]["is_contaminant"], 0)

    def test_calls_are_driven_by_prevalence_not_abundance(self):
        """An abundant feature absent from the blanks must never be flagged."""
        from databases.database_metrics import decontam_prevalence

        with tempfile.TemporaryDirectory() as d:
            matrix = _write(
                Path(d) / "m.csv",
                "#OTU ID,c1,c2,c3,s1,s2,s3\n"
                "tiny_but_in_blanks,1,1,1,0,0,0\n"
                "huge_but_clean,0,0,0,100000,90000,80000\n",
            )
            rows = {row["feature"]: row for row in decontam_prevalence(matrix, ["c1", "c2", "c3"])}

        self.assertGreater(rows["huge_but_clean"]["score"], rows["tiny_but_in_blanks"]["score"])

    def test_controls_only_run_yields_nothing(self):
        from databases.database_metrics import decontam_prevalence

        with tempfile.TemporaryDirectory() as d:
            matrix = _write(Path(d) / "m.csv", "#OTU ID,c1,c2\nf1,1,1\n")
            self.assertEqual(decontam_prevalence(matrix, ["c1", "c2"]), [])


class TaxonPrevalenceTests(unittest.TestCase):
    def test_controls_are_excluded_from_the_denominator(self):
        from databases.database_metrics import taxon_prevalence

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            rows = {
                row["taxon"]: row
                for row in taxon_prevalence(_matrix(root / "m.csv"), _consensus(root / "c.csv"), exclude_samples=["ctrl1", "ctrl2"])
            }

        self.assertEqual(rows["Diplostomum"]["n_samples"], 4)
        self.assertEqual(rows["Diplostomum"]["n_detected"], 4)
        self.assertEqual(rows["Diplostomum"]["prevalence"], 1.0)
        # Even at 4/4 the interval acknowledges four fish is not proof.
        self.assertLess(rows["Diplostomum"]["ci_low"], 1.0)

    def test_features_without_a_taxon_collapse_into_unassigned(self):
        from databases.database_metrics import taxon_prevalence

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            rows = {row["taxon"] for row in taxon_prevalence(_matrix(root / "m.csv"), _consensus(root / "c.csv"))}

        self.assertIn("Unassigned", rows)


class PermanovaTests(unittest.TestCase):
    def test_two_groups_of_two_are_testable(self):
        from databases.database_metrics import permanova_by_group

        with tempfile.TemporaryDirectory() as d:
            result = permanova_by_group(_matrix(Path(d) / "m.csv"), {"s1": "A", "s2": "A", "s3": "B", "s4": "B"})

        self.assertEqual(result["n_groups"], 2)
        self.assertEqual(result["n_samples"], 4)
        self.assertIn("p_value", result)

    def test_a_singleton_group_tests_nothing(self):
        from databases.database_metrics import permanova_by_group

        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(permanova_by_group(_matrix(Path(d) / "m.csv"), {"s1": "A", "s2": "B", "s3": "C"}), {})

    def test_no_groups_returns_nothing(self):
        from databases.database_metrics import permanova_by_group

        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(permanova_by_group(_matrix(Path(d) / "m.csv"), {}), {})


class BlastCoverageTests(unittest.TestCase):
    def test_ncbi_outfmt_requests_query_coverage(self):
        from databases.ncbi_db import NCBIConfig

        self.assertIn("qcovhsp", NCBIConfig().outfmt)

    def test_custom_blast_besthit_carries_coverage(self):
        from databases.custom_db import _BESTHIT_COLUMNS

        self.assertIn("qcovhsp", _BESTHIT_COLUMNS)

    def test_older_raw_tsv_without_qcovhsp_still_parses(self):
        from databases.custom_db import _parse_best_hits

        with tempfile.TemporaryDirectory() as d:
            tsv = _write(Path(d) / "raw.tsv", "q1\ts1\t99.0\t300\t0\t0\t1\t300\t1\t300\t1e-50\t500\n")
            best = _parse_best_hits(tsv)

        self.assertEqual(best["q1"]["pident"], "99.0")
        self.assertEqual(best["q1"]["qcovhsp"], "")


if __name__ == "__main__":
    unittest.main()


class DatabaseFigureTests(unittest.TestCase):
    def _frame(self, rows):
        import pandas as pd

        return pd.DataFrame(rows)

    def test_confidence_panels_omit_evidence_the_sources_never_produce(self):
        """A SINTAX-only run must not draw an empty identity axis."""
        from ui.plots import confidence_histogram_figure

        sintax = self._frame([{"identity": "", "coverage": "", "confidence": 0.9}])
        titles = [annotation.text for annotation in confidence_histogram_figure(sintax).layout.annotations]
        self.assertEqual(titles, ["Confidence / bootstrap"])

    def test_confidence_panels_appear_when_an_aligning_source_ran(self):
        from ui.plots import confidence_histogram_figure

        mixed = self._frame([{"identity": 99.0, "coverage": 98.0, "confidence": 820.0}, {"identity": "", "coverage": "", "confidence": 0.9}])
        titles = [annotation.text for annotation in confidence_histogram_figure(mixed).layout.annotations]
        self.assertEqual(titles, ["% identity", "% query coverage", "Confidence / bootstrap"])

    def test_sunburst_parent_equals_the_sum_of_its_children(self):
        from ui.plots import taxonomy_sunburst_figure

        consensus = self._frame(
            [
                {"Consensus_Kingdom": "Eukaryota", "Consensus_Phylum": "Platyhelminthes", "Consensus_Class": "", "Consensus_Order": "", "Consensus_Family": "", "Consensus_Genus": ""},
                {"Consensus_Kingdom": "Eukaryota", "Consensus_Phylum": "", "Consensus_Class": "", "Consensus_Order": "", "Consensus_Family": "", "Consensus_Genus": ""},
            ]
        )
        trace = taxonomy_sunburst_figure(consensus).data[0]
        roots = [value for node, value in zip(trace.ids, trace.values) if "/" not in node]
        self.assertEqual(roots, [2])
        self.assertEqual(trace.branchvalues, "total")

    def test_prevalence_error_bars_are_asymmetric_around_the_estimate(self):
        from ui.plots import prevalence_figure

        figure = prevalence_figure(self._frame([{"taxon": "D", "prevalence": 0.6, "ci_low": 0.2306, "ci_high": 0.8824, "total_reads": 10}]))
        self.assertAlmostEqual(float(figure.data[0].error_x.array[0]), 28.24, places=1)
        self.assertAlmostEqual(float(figure.data[0].error_x.arrayminus[0]), 36.94, places=1)

    def test_empty_frames_render_nothing_rather_than_an_empty_axis(self):
        from ui.plots import contaminant_figure, prevalence_figure, taxonomic_resolution_figure

        empty = self._frame([])
        self.assertIsNone(taxonomic_resolution_figure(empty))
        self.assertIsNone(prevalence_figure(empty))
        self.assertIsNone(contaminant_figure(empty))


class PhyloseqBundleTests(unittest.TestCase):
    """The phyloseq export exists because BIOM drops the sample metadata.

    phyloseq intersects its three tables on their ids, so a mismatch does not raise —
    it silently yields an object with fewer taxa or fewer samples than the run had.
    These tests pin the alignment and the orientation rather than just the file list.
    """

    def _bundle(self, root: Path, **kwargs):
        from databases.database_metrics import write_phyloseq_bundle

        return write_phyloseq_bundle(
            _matrix(root / "m.csv"),
            _consensus(root / "c.csv"),
            root / "phyloseq",
            **kwargs,
        )

    def _rows(self, path: Path) -> list[list[str]]:
        # Split on newlines only: .strip() would eat the trailing tabs of the last row,
        # which are exactly the empty ranks these tests are checking for.
        return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines() if line]

    def test_the_bundle_is_the_three_tables_plus_the_script(self):
        with tempfile.TemporaryDirectory() as d:
            written = self._bundle(Path(d))
            self.assertEqual(
                sorted(written),
                ["phyloseq_otu_table", "phyloseq_sample_data", "phyloseq_script", "phyloseq_tax_table"],
            )
            self.assertTrue(all(path.exists() for path in written.values()))

    def test_otu_table_is_features_by_samples(self):
        # taxa_are_rows = TRUE in the script: transposing this silently swaps what a
        # sample means for every downstream call.
        with tempfile.TemporaryDirectory() as d:
            rows = self._rows(self._bundle(Path(d))["phyloseq_otu_table"])
            self.assertEqual(rows[0], ["feature_id", "ctrl1", "ctrl2", "s1", "s2", "s3", "s4"])
            self.assertEqual([r[0] for r in rows[1:]], ["to_species", "to_genus", "to_kingdom"])
            self.assertEqual(rows[1], ["to_species", "0", "0", "100", "120", "90", "80"])

    def test_feature_ids_line_up_across_otu_and_tax_tables(self):
        with tempfile.TemporaryDirectory() as d:
            written = self._bundle(Path(d))
            otu_ids = [r[0] for r in self._rows(written["phyloseq_otu_table"])[1:]]
            tax_ids = [r[0] for r in self._rows(written["phyloseq_tax_table"])[1:]]
            self.assertEqual(otu_ids, tax_ids)

    def test_sample_ids_line_up_across_otu_and_sample_data(self):
        with tempfile.TemporaryDirectory() as d:
            written = self._bundle(Path(d))
            otu_samples = self._rows(written["phyloseq_otu_table"])[0][1:]
            sample_ids = [r[0] for r in self._rows(written["phyloseq_sample_data"])[1:]]
            self.assertEqual(otu_samples, sample_ids)

    def test_an_unassigned_rank_is_blank_not_a_fabricated_taxon(self):
        # "unclassified" here would be treated by tax_glom() as one real taxon and
        # merge unrelated features under it.
        with tempfile.TemporaryDirectory() as d:
            rows = {r[0]: r for r in self._rows(self._bundle(Path(d))["phyloseq_tax_table"])[1:]}
            self.assertEqual(rows["to_species"][7], "Diplostomum spathaceum")
            self.assertEqual(rows["to_genus"][7], "")  # genus known, species not
            self.assertEqual(rows["to_genus"][6], "Ralstonia")
            self.assertEqual(rows["to_kingdom"][2:], [""] * 6)

    def test_groups_and_controls_survive_into_sample_data(self):
        # This is the whole reason the bundle exists: BIOM writes columns metadata as
        # null, so exporting that way loses both of these.
        with tempfile.TemporaryDirectory() as d:
            written = self._bundle(
                Path(d),
                control_samples=["ctrl1", "ctrl2"],
                groups={"s1": "infected", "s2": "infected"},
            )
            rows = {r[0]: r for r in self._rows(written["phyloseq_sample_data"])[1:]}
            self.assertEqual(rows["ctrl1"][2], "yes")
            self.assertEqual(rows["s1"][1:3], ["infected", "no"])
            # A sample with no group is labelled, not left empty: phyloseq drops NA
            # levels from a factor, which would quietly remove those samples from plots.
            self.assertEqual(rows["s3"][1], "ungrouped")

    def test_read_totals_are_the_column_sums(self):
        with tempfile.TemporaryDirectory() as d:
            rows = {r[0]: r for r in self._rows(self._bundle(Path(d))["phyloseq_sample_data"])[1:]}
            self.assertEqual(rows["s1"][3], "108")  # 100 + 5 + 3
            self.assertEqual(rows["ctrl1"][3], "10")

    def test_counts_are_raw_not_relative(self):
        # Rarefying/normalising is the analyst's call; phyloseq has functions for it.
        with tempfile.TemporaryDirectory() as d:
            rows = self._rows(self._bundle(Path(d))["phyloseq_otu_table"])
            self.assertEqual(rows[1][3], "100")

    def test_a_missing_matrix_yields_no_bundle(self):
        with tempfile.TemporaryDirectory() as d:
            from databases.database_metrics import write_phyloseq_bundle

            self.assertEqual(write_phyloseq_bundle(Path(d) / "absent.csv", None, Path(d) / "ps"), {})

    def test_a_run_without_taxonomy_still_exports_counts_and_metadata(self):
        # No consensus is not a reason to withhold the OTU table; the tax table is
        # simply empty and phyloseq still loads.
        with tempfile.TemporaryDirectory() as d:
            from databases.database_metrics import write_phyloseq_bundle

            written = write_phyloseq_bundle(_matrix(Path(d) / "m.csv"), None, Path(d) / "ps")
            rows = self._rows(written["phyloseq_tax_table"])
            self.assertEqual(len(rows), 4)
            self.assertEqual(rows[1][1:], [""] * 7)
