"""Fixes found while running the Guadiana dataset (F5).

Each test fails on commit 7765c4f and passes after the fix. None of the fixes changes a
feature table, a taxonomic assignment or a tree: they repair a crash (UniFrac), a wrong
status label (custom BLAST), dropped tips (Faith's PD), a silent no-op (COI on "auto")
and parameter keys that happened to match their defaults (FastTree).
"""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

from core.models import AnalysisMode, Marker, Sample, SamplePaths


class UnifracOrdinationTests(unittest.TestCase):
    """scikit-bio >= 0.7 rejects ``otu_ids=``: UniFrac raised "taxa is required"."""

    def test_weighted_unifrac_ordinates_with_current_scikit_bio(self):
        from analysis.ordination import ordinate

        rng = np.random.default_rng(1)
        a = rng.multinomial(2_000, [0.4, 0.3, 0.2, 0.05, 0.05], size=3)
        b = rng.multinomial(2_000, [0.05, 0.05, 0.2, 0.3, 0.4], size=3)
        counts = np.vstack([a, b]).astype(float)
        samples = ["A0", "A1", "A2", "B0", "B1", "B2"]
        features = [f"OTU{i}" for i in range(1, 6)]
        with tempfile.TemporaryDirectory() as d:
            tree = Path(d) / "tree.nwk"
            tree.write_text("(((OTU1:0.1,OTU2:0.2):0.3,OTU3:0.4):0.2,(OTU4:0.1,OTU5:0.2):0.5);\n")
            result = ordinate(counts, samples, features, metric="weighted_unifrac", tree_path=tree)
        self.assertEqual(len(result["coordinates"] if "coordinates" in result else result["samples"]), 6)


class CustomBlastStatusTests(unittest.TestCase):
    """A custom BLAST with matches was reported as ``no_hits`` on every run."""

    def _reconcile(self, matched):
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "custom_blast_raw.tsv"
            raw.write_text("OTU1\tsub\t99\t300\t0\t0\t1\t300\t1\t300\t1e-50\t500\t100\n")
            log = Path(d) / "custom_blast.log"
            log.write_text("Custom BLAST: 1/2 OTUs matched.\n")
            entry = {"source": "Custom database (BLAST)", "status": "ok", "detail": "x",
                     "key": "CUSTOM", "raw": raw, "log": log, "matched": matched}
            return _reconcile_source_status([entry], {}, 2, marker="COI")[0]

    def test_matches_reported_by_blast_count_as_assignments(self):
        out = self._reconcile(1)
        self.assertEqual(out["status"], "ok")
        self.assertTrue(out["detail"].startswith("1/2"))

    def test_zero_matches_is_still_no_hits(self):
        self.assertEqual(self._reconcile(0)["status"], "no_hits")

    def test_the_hint_never_reaches_the_csv(self):
        self.assertNotIn("matched", self._reconcile(1))


class FaithReversedTipTests(unittest.TestCase):
    """MAFFT --adjustdirection renames reversed records ``_R_<id>``; those tips were lost."""

    def test_reversed_tip_contributes_its_branch(self):
        from phylo.phylo_diversity import faith_pd

        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            tree = d / "tree.nwk"
            tree.write_text("((ASV1:1.0,_R_ASV2:2.0):0.5,ASV3:1.0);\n")
            mapping = d / "reversed_asvs.csv"
            with mapping.open("w", newline="") as handle:
                w = csv.writer(handle)
                w.writerow(["ASV_ID", "Original_Header"])
                w.writerows([["ASV1", "OTU1"], ["ASV2", "OTU2"], ["ASV3", "OTU3"]])
            matrix = d / "otu_matrix.csv"
            matrix.write_text("#OTU ID,S1\nOTU1,0\nOTU2,5\nOTU3,0\n")
            result = faith_pd(tree, matrix, mapping)
        self.assertEqual(result["matched_features"], 1)
        self.assertGreater(result["pd"][0], 0.0)


class AutoPrimerWarningTests(unittest.TestCase):
    """A COI run left on "auto" removed no primers and said nothing."""

    def _config(self, marker, params):
        sample = Sample("s1", "S1", SamplePaths(r1=Path("r1.fq.gz"), r2=None), AnalysisMode.SE, marker)
        return SimpleNamespace(samples=[sample], vsearch_params=params, output_dir=Path("."))

    def test_coi_on_auto_is_named(self):
        from preprocess.cutadapt_runner import auto_without_pair

        self.assertEqual(auto_without_pair(self._config(Marker.COI, {"primer_preset": "auto"})), "COI")

    def test_18s_on_auto_has_a_pair(self):
        from preprocess.cutadapt_runner import auto_without_pair

        self.assertEqual(auto_without_pair(self._config(Marker.SSU_18S, {"primer_preset": "auto"})), "")

    def test_explicit_opt_out_is_not_a_warning(self):
        from preprocess.cutadapt_runner import auto_without_pair

        self.assertEqual(auto_without_pair(self._config(Marker.COI, {"primer_preset": ""})), "")

    def test_passthrough_carries_the_warning(self):
        from preprocess import cutadapt_runner as cr

        out = cr.run_cutadapt(self._config(Marker.COI, {"primer_preset": "auto"}).samples[0],
                              self._config(Marker.COI, {"primer_preset": "auto"}))
        self.assertTrue(out.success)
        self.assertIn("COI", out.error_message)
        self.assertEqual(out.metrics["primers_removed"], 0.0)


class FastTreeKeyTests(unittest.TestCase):
    """RunConfig stores ``fasttree_gtr``/``fasttree_gamma``; the runner read ``gtr``/``gamma``."""

    def test_runconfig_keys_reach_the_command(self):
        from phylo import fasttree_runner as ft

        with tempfile.TemporaryDirectory() as d:
            aligned = Path(d) / "aligned.fasta"
            aligned.write_text("".join(f">S{i}\nACGTACGTAC\n" for i in range(5)))
            seen = {}

            def fake_run(command, **kwargs):
                seen["command"] = command
                (Path(d) / "out" / "tree.nwk").write_text("(S0,S1,(S2,(S3,S4)));\n")
                return SimpleNamespace(returncode=0, stdout="(S0,S1,(S2,(S3,S4)));\n", stderr="")

            (Path(d) / "out").mkdir()
            with mock.patch.object(ft.subprocess, "run", side_effect=fake_run):
                ft.run_fasttree(aligned, Path(d) / "out", {"fasttree_gtr": True, "fasttree_gamma": False})
        self.assertIn("-gtr", seen["command"])
        self.assertNotIn("-gamma", seen["command"])


if __name__ == "__main__":
    unittest.main()
