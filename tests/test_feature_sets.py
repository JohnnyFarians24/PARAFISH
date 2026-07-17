"""OTU and ASV feature sets: generated alike, reported apart."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from analysis.analysis_metrics import build_metrics_from_outputs, feature_occupancy, feature_sets_from_outputs
from analysis.feature_set_comparison import build_summary_rows, build_taxon_bridge
from databases.taxonomy_runner import (
    _build_consensus_row,
    _build_hits_rows,
    _genus_from_binomial,
    best_alignment_record,
)

import numpy as np


class FeatureSetDetectionTests(unittest.TestCase):
    def _outputs(self, root: Path, *, otu: bool, asv: bool) -> dict:
        outputs: dict = {}
        if otu:
            fasta = root / "otus.fasta"
            fasta.write_text(">OTU1;size=5\nACGT\n", encoding="utf-8")
            matrix = root / "otu_matrix.csv"
            matrix.write_text("OTU_ID,S1,S2\nOTU1,3,2\n", encoding="utf-8")
            outputs["combined_fasta_vsearch"] = str(fasta)
            outputs["otu_matrix"] = str(matrix)
        if asv:
            fasta = root / "asv.fasta"
            fasta.write_text(">ASV1;size=4\nACGT\n", encoding="utf-8")
            matrix = root / "asv_abundance.csv"
            matrix.write_text("ASV_ID,S1,S2\nASV1,4,0\n", encoding="utf-8")
            outputs["combined_fasta_dada2"] = str(fasta)
            outputs["dada2"] = {"files": {"asv_sequences_fasta": str(fasta), "asv_abundance_csv": str(matrix)}}
        return outputs

    def test_both_branches_yield_two_independent_feature_sets(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            sets = feature_sets_from_outputs(self._outputs(root, otu=True, asv=True))
        self.assertEqual([s["key"] for s in sets], ["otu", "asv"])

    def test_metrics_are_written_to_separate_directories(self):
        """The whole point: an OTU table and an ASV table can never be the same file."""
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            written = build_metrics_from_outputs(root, self._outputs(root, otu=True, asv=True))

            self.assertIn("otu.prevalence_csv", written)
            self.assertIn("asv.prevalence_csv", written)
            self.assertNotEqual(written["otu.prevalence_csv"], written["asv.prevalence_csv"])
            self.assertTrue(Path(written["otu.prevalence_csv"]).is_relative_to(root / "metrics" / "otu"))
            self.assertTrue(Path(written["asv.prevalence_csv"]).is_relative_to(root / "metrics" / "asv"))

    def test_dada2_only_run_still_produces_a_feature_set(self):
        """DADA2 alone must reach the Databases stage; it used to be dropped entirely."""
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            sets = feature_sets_from_outputs(self._outputs(root, otu=False, asv=True))
        self.assertEqual([s["key"] for s in sets], ["asv"])


class ComparisonTests(unittest.TestCase):
    def test_summary_reports_one_row_per_set_never_a_merged_row(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            otu_matrix = root / "otu.csv"
            otu_matrix.write_text("OTU_ID,S1,S2\nOTU1,3,2\n", encoding="utf-8")
            asv_matrix = root / "asv.csv"
            asv_matrix.write_text("ASV_ID,S1,S2\nASV1,4,0\nASV2,1,1\n", encoding="utf-8")
            rows = build_summary_rows(
                [
                    {"key": "otu", "label": "OTU", "matrix": otu_matrix},
                    {"key": "asv", "label": "ASV", "matrix": asv_matrix},
                ]
            )
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["n_features"], 1)
        self.assertEqual(rows[1]["n_features"], 2)
        # Counts stay on their own row; nothing is summed across sets.
        self.assertEqual({row["feature_set"] for row in rows}, {"otu", "asv"})

    def test_taxon_bridge_joins_on_taxon_and_keeps_counts_apart(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            otu_c = root / "otu_consensus.csv"
            otu_c.write_text("ASV_ID,Consensus_Genus\nOTU1,Diplostomum\n", encoding="utf-8")
            asv_c = root / "asv_consensus.csv"
            asv_c.write_text("ASV_ID,Consensus_Genus\nASV1,Diplostomum\nASV2,Anisakis\n", encoding="utf-8")
            rows = build_taxon_bridge(otu_c, asv_c, None, None, rank="genus")

        by_taxon = {row["taxon"]: row for row in rows}
        shared = by_taxon["Diplostomum"]
        self.assertEqual(shared["agreement"], "both")
        # One OTU and one ASV both found it — reported separately, not as "2 features".
        self.assertEqual((shared["n_otus"], shared["n_asvs"]), (1, 1))
        self.assertEqual(by_taxon["Anisakis"]["agreement"], "asv_only")


class OccupancyTests(unittest.TestCase):
    def test_reports_which_samples_not_just_how_many(self):
        counts = np.array([[5, 0, 2], [3, 0, 0]], dtype=float)  # samples x features
        rows = feature_occupancy(counts, ["F1", "F2", "F3"], ["S1", "S2"])
        by_feature = {row["feature"]: row for row in rows}

        self.assertEqual(by_feature["F1"]["samples_present"], "S1;S2")
        self.assertEqual(by_feature["F1"]["n_samples_present"], 2)
        self.assertEqual(by_feature["F1"]["is_shared"], 1)
        self.assertEqual(by_feature["F3"]["samples_present"], "S1")
        self.assertEqual(by_feature["F3"]["is_unique"], 1)
        self.assertNotIn("F2", by_feature)  # absent everywhere
        # Shared features come first: they are the ones a cross-sample study is about.
        self.assertEqual(rows[0]["feature"], "F1")


class EvidenceColumnTests(unittest.TestCase):
    def test_consensus_carries_identity_coverage_and_reliability(self):
        """These used to exist only in source_assignments.csv, which nobody opened."""
        records = [
            {"source": "NCBI", "species": "Diplostomum spathaceum", "genus": "Diplostomum", "similarity": "99.2", "coverage": "98"},
            {"source": "SILVA", "genus": "Diplostomum", "confidence": "0.9"},
        ]
        row = _build_consensus_row("OTU1", records, "lca", "OTU")

        self.assertEqual(row["feature_type"], "OTU")
        self.assertEqual(row["Best_Identity"], "99.2")
        self.assertEqual(row["Best_Coverage"], "98")
        self.assertEqual(row["Best_Source"], "NCBI")
        self.assertEqual(row["Reliability"], "high")

    def test_identity_below_the_rank_threshold_is_not_called_reliable(self):
        records = [{"source": "NCBI", "species": "Diplostomum spathaceum", "similarity": "85.0", "coverage": "99"}]
        row = _build_consensus_row("OTU1", records, "lca", "OTU")
        # 85% cannot support a species-level call (threshold is 97%).
        self.assertEqual(row["Reliability"], "low")

    def test_best_record_prefers_identity_then_coverage(self):
        records = [
            {"source": "NCBI", "species": "Diplostomum spathaceum", "similarity": "97.0", "coverage": "60"},
            {"source": "BOLD", "species": "Diplostomum pseudospathaceum", "similarity": "99.0", "coverage": "50"},
            {"source": "SILVA", "species": "Diplostomum sp.", "confidence": "0.99"},  # no identity: not comparable
        ]
        best = best_alignment_record(records)
        self.assertEqual(best["source"], "BOLD")

    def test_a_source_that_named_nothing_is_not_the_best_hit(self):
        """BOLD writes its non-identifications as rows with similarity=0.0 — not None, so
        they used to score. With SILVA/PR2 excluded for being SINTAX, BOLD's non-answer was
        the only scoring record, and every feature was credited Best_Source=BOLD,
        Best_Identity=0.0 while the species on the row had come from SILVA."""
        records = [
            {"source": "BOLD", "species": "", "similarity": "0.0", "flags": "||||"},
            {"source": "SILVA", "species": "Gadus morhua", "confidence": "1.00"},
        ]
        self.assertIsNone(best_alignment_record(records))

        row = _build_consensus_row("OTU1", records, "lca", "OTU")
        self.assertNotEqual(row["Best_Source"], "BOLD")
        self.assertEqual(row["Best_Source"], "SILVA")  # falls through to the SINTAX evidence

    def test_a_zero_identity_alignment_is_not_a_hit(self):
        records = [{"source": "NCBI", "species": "Something", "similarity": "0.0", "coverage": "0"}]
        self.assertIsNone(best_alignment_record(records))

    def test_hits_are_ranked_best_first_so_the_runner_up_stays_reachable(self):
        hits = {
            "OTU1": [
                {"source": "NCBI", "species": "Uncultured", "similarity": "91.0", "coverage": "80"},
                {"source": "NCBI", "species": "Diplostomum spathaceum", "similarity": "99.0", "coverage": "95"},
            ]
        }
        rows = _build_hits_rows(hits, "OTU")
        self.assertEqual([row["hit_rank"] for row in rows], [1, 2])
        self.assertEqual(rows[0]["identity"], "99.0")
        self.assertEqual(rows[1]["identity"], "91.0")  # kept, not discarded
        self.assertTrue(all(row["feature_type"] == "OTU" for row in rows))


class GenusRecoveryTests(unittest.TestCase):
    def test_genus_is_recovered_from_a_blast_binomial(self):
        """BLAST returns only a scientific name, so the genus column arrived empty."""
        self.assertEqual(_genus_from_binomial("Micromesistius poutassou"), "Micromesistius")
        self.assertEqual(_genus_from_binomial("Diplostomum sp."), "Diplostomum")

    def test_placeholder_names_do_not_invent_a_genus(self):
        for name in ("uncultured eukaryote", "Eukaryota environmental sample", "Bacteria", "unidentified marine organism"):
            self.assertEqual(_genus_from_binomial(name), "", name)


if __name__ == "__main__":
    unittest.main()
