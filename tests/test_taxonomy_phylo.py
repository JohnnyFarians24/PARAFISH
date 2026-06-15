"""Taxonomy and phylogenetic contract tests."""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from core.models import AnalysisBranch, AnalysisMode, Marker, PhyloMode, RunConfig, Sample, SamplePaths, StepOutput
from databases.taxonomy_runner import build_consensus, run_sintax
from phylo.phylo_runner import run_phylo


class TaxonomyPhyloTests(unittest.TestCase):
    def _sample(self, temp_dir: Path) -> Sample:
        r1 = temp_dir / "sample_R1.fastq"
        r1.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
        return Sample(
            internal_id="sample_001",
            display_label="Sample 1",
            paths=SamplePaths(r1=r1),
            mode=AnalysisMode.SE,
            marker=Marker.COI,
        )

    def test_run_sintax_output_exists(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta = temp_path / "asv_sequences.fasta"
            fasta.write_text(">ASV1\n" + "A" * 120 + "\n", encoding="utf-8")
            db = temp_path / "db.fasta"
            db.write_text(">db\n" + "A" * 120 + "\n", encoding="utf-8")

            def fake_run(command, capture_output, text, check):
                del capture_output, text, check
                result_index = command.index("--tabbedout") + 1
                Path(command[result_index]).write_text("ASV1\tTaxon\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, stdout="sintax ok", stderr="")

            with patch("databases.taxonomy_runner.subprocess.run", side_effect=fake_run):
                result = run_sintax(fasta, db, 0.8, 8)

            self.assertTrue(result.exists())

    def test_build_consensus_merges_rows(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bold = temp_path / "bold_result.parquet"
            bold.write_text("ASV_ID,BOLD_species,BOLD_similarity\nASV1,SpeciesA,0.97\n", encoding="utf-8")
            sintax = temp_path / "sintax_result.csv"
            sintax.write_text("ASV_ID,SINTAX_species,SINTAX_bootstrap\nASV1,SpeciesA,0.91\n", encoding="utf-8")

            table = build_consensus(bold, sintax, None, method="lca")
            self.assertEqual(len(table), 1)
            self.assertEqual(table.rows[0]["Source_of_consensus"], "lca")

    def test_run_phylo_fast_routes_to_fasttree(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            aligned = temp_path / "aligned.fasta"
            aligned.write_text(">ASV1\n" + "A" * 120 + "\n>ASV2\n" + "T" * 120 + "\n>ASV3\n" + "G" * 120 + "\n>ASV4\n" + "C" * 120 + "\n", encoding="utf-8")
            sample = self._sample(temp_path)
            config = RunConfig(run_id="run_phylo_fast", samples=[sample], branch=AnalysisBranch.VSEARCH, phylo_mode=PhyloMode.FAST, output_dir=temp_path / "out")

            mafft_result = StepOutput(step_name="mafft", success=True, files={"aligned_fasta": aligned, "log": temp_path / "mafft.log", "reversed_asvs_csv": temp_path / "rev.csv"})
            fasttree_result = StepOutput(step_name="fasttree", success=True, files={"tree_nwk": temp_path / "tree.nwk", "log": temp_path / "fasttree.log"})

            with patch("phylo.phylo_runner.run_mafft", return_value=mafft_result) as mocked_mafft, patch("phylo.phylo_runner.run_fasttree", return_value=fasttree_result) as mocked_fasttree, patch("phylo.phylo_runner.run_raxml") as mocked_raxml:
                result = run_phylo(aligned, config, temp_path / "phylo")

            self.assertTrue(result.success)
            mocked_mafft.assert_called_once()
            mocked_fasttree.assert_called_once()
            mocked_raxml.assert_not_called()

    def test_run_phylo_hybrid_routes_to_raxml(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            aligned = temp_path / "aligned.fasta"
            aligned.write_text(">ASV1\n" + "A" * 120 + "\n>ASV2\n" + "T" * 120 + "\n>ASV3\n" + "G" * 120 + "\n>ASV4\n" + "C" * 120 + "\n", encoding="utf-8")
            sample = self._sample(temp_path)
            config = RunConfig(run_id="run_phylo_hybrid", samples=[sample], branch=AnalysisBranch.VSEARCH, phylo_mode=PhyloMode.HYBRID, output_dir=temp_path / "out")

            mafft_result = StepOutput(step_name="mafft", success=True, files={"aligned_fasta": aligned, "log": temp_path / "mafft.log", "reversed_asvs_csv": temp_path / "rev.csv"})
            tree_nwk = temp_path / "tree.nwk"
            tree_nwk.write_text("(ASV1:0.1,ASV2:0.2,ASV3:0.3,ASV4:0.4);", encoding="utf-8")
            fasttree_result = StepOutput(step_name="fasttree", success=True, files={"tree_nwk": tree_nwk, "log": temp_path / "fasttree.log"})
            raxml_result = StepOutput(step_name="raxml", success=True, files={"best_tree": temp_path / "bestTree", "support_tree": temp_path / "support", "bootstraps": temp_path / "boots", "best_model": temp_path / "bestModel", "log": temp_path / "raxml.log"})

            with patch("phylo.phylo_runner.run_mafft", return_value=mafft_result), patch("phylo.phylo_runner.run_fasttree", return_value=fasttree_result), patch("phylo.phylo_runner.run_raxml", return_value=raxml_result) as mocked_raxml:
                result = run_phylo(aligned, config, temp_path / "phylo")

            self.assertTrue(result.success)
            mocked_raxml.assert_called_once()


if __name__ == "__main__":
    unittest.main()