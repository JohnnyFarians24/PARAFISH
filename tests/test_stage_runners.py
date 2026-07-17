"""Runner contract tests for preprocessing, analysis, taxonomy and phylogenetics."""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from analysis.dada2_runner import run_dada2
from analysis.vsearch_runner import run_vsearch_filter_sample, run_vsearch_pool_and_cluster
from core.utils import taxonomy_fields_to_sintax
from core.models import AnalysisBranch, AnalysisMode, Marker, PhyloMode, RunConfig, Sample, SamplePaths, StepOutput
from core.validation import validate_fasta_headers
from databases.taxonomy_runner import run_boldigger
from phylo.phylo_runner import run_fasttree, run_mafft, run_raxml, sanitize_fasta_for_phylo
from preprocess.pear_runner import run_pear
from preprocess.trimmomatic_runner import run_trimmomatic


class RunnerContractTests(unittest.TestCase):
    def _write_fastq(self, path: Path, reads: int = 1) -> Path:
        lines = []
        for index in range(reads):
            lines.extend([f"@read{index + 1}", "ACGT", "+", "!!!!"])
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def _sample(self, temp_dir: Path, sample_id: str = "sample_001", mode: AnalysisMode = AnalysisMode.PE, marker: Marker = Marker.COI) -> Sample:
        r1 = self._write_fastq(temp_dir / f"{sample_id}_R1.fastq", reads=2)
        r2 = self._write_fastq(temp_dir / f"{sample_id}_R2.fastq", reads=2) if mode == AnalysisMode.PE else None
        return Sample(
            internal_id=sample_id,
            display_label="Sample 1",
            paths=SamplePaths(r1=r1, r2=r2),
            mode=mode,
            marker=marker,
        )

    def test_pear_rejects_dada2_branch(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            sample = self._sample(temp_path)
            config = RunConfig(run_id="run1", samples=[sample], branch=AnalysisBranch.DADA2, output_dir=temp_path / "out")
            trimmed = StepOutput(step_name="trimmomatic", success=True, files={"r1_trimmed": sample.paths.r1, "r2_trimmed": sample.paths.r2})

            with self.assertRaises(RuntimeError):
                run_pear(sample, config, trimmed)

    def test_trimmomatic_pe_success(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            sample = self._sample(temp_path)
            config = RunConfig(run_id="run2", samples=[sample], branch=AnalysisBranch.VSEARCH, output_dir=temp_path / "out")

            def fake_run(command, capture_output, text, check):
                del capture_output, text, check
                for arg in command:
                    if str(arg).endswith("_R1_trimmed.fastq.gz"):
                        self._write_fastq(Path(arg), reads=2)
                    if str(arg).endswith("_R2_trimmed.fastq.gz"):
                        self._write_fastq(Path(arg), reads=2)
                return subprocess.CompletedProcess(command, 0, stdout="trim ok", stderr="")

            with patch("preprocess.trimmomatic_runner.subprocess.run", side_effect=fake_run):
                result = run_trimmomatic(sample, config)

            self.assertTrue(result.success)
            self.assertIn("r1_trimmed", result.files)

    def test_vsearch_filter_sample_tags_reads(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            sample = self._sample(temp_path)
            config = RunConfig(run_id="run3a", samples=[sample], branch=AnalysisBranch.VSEARCH, output_dir=temp_path / "out")
            assembled_path = temp_path / "assembled.fastq"
            assembled_path.write_text("@r1\n" + "A" * 220 + "\n+\n" + ("!" * 220) + "\n", encoding="utf-8")

            def fake_run(command, capture_output, text, check):
                del capture_output, text, check
                joined = " ".join(str(arg) for arg in command)
                if "--fastq_filter" in joined:
                    out_index = command.index("--fastaout") + 1
                    Path(command[out_index]).write_text(">r1;sample=sample_001;\n" + "A" * 220 + "\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

            with patch("analysis.vsearch_runner.subprocess.run", side_effect=fake_run):
                result = run_vsearch_filter_sample(assembled_path, config, temp_path / "out" / "sample_001", "sample_001")

            self.assertTrue(result.success)
            self.assertIn("filtered_fasta", result.files)

    def test_vsearch_pool_fails_empty_otu(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            config = RunConfig(run_id="run3b", samples=[], branch=AnalysisBranch.VSEARCH, output_dir=temp_path / "out")
            filtered = temp_path / "sample_001_filtered.fasta"
            filtered.write_text(">r1;sample=sample_001;\n" + "A" * 220 + "\n", encoding="utf-8")

            def fake_run(command, capture_output, text, check):
                del capture_output, text, check
                return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

            with patch("analysis.vsearch_runner.subprocess.run", side_effect=fake_run):
                result = run_vsearch_pool_and_cluster({"sample_001": filtered}, config, temp_path / "pool_out")

            self.assertFalse(result.success)
            self.assertIn("0 OTUs", result.error_message)

    def test_dada2_runs_without_refdb(self):
        """No reference is needed to call ASVs — taxonomy comes from the shared stage.

        Requiring ``ref_fasta`` up front is what made every DADA2 run die with
        "DADA2 requires a reference FASTA database" before denoising a single read,
        even though the R script has always handled its absence.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            sample = self._sample(temp_path)
            config = RunConfig(run_id="run4", samples=[sample], branch=AnalysisBranch.DADA2, output_dir=temp_path / "out")
            self.assertIsNone(config.dada2_params.get("ref_fasta"))
            # Must not raise. Whether R is installed decides success/failure here; the
            # contract under test is only that the missing reference is not fatal.
            result = run_dada2([sample], config)
            self.assertNotIn("requires a reference", str(result.error_message or ""))

    @patch("core.preflight.dada2_memory_shortfall", lambda *a, **k: None)
    def test_dada2_success_creates_outputs(self):
        # The memory gate is asserted in test_partial_analysis; here it must not turn a
        # contract test into something that passes or fails with the machine's free RAM.
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            sample = self._sample(temp_path, marker=Marker.SSU_18S)
            config = RunConfig(run_id="run5", samples=[sample], branch=AnalysisBranch.DADA2, output_dir=temp_path / "out")
            ref_fasta = temp_path / "ref.fasta"
            ref_fasta.write_text(">ref\nACGT\n", encoding="utf-8")
            config.dada2_params["ref_fasta"] = ref_fasta

            def fake_run(command, capture_output, text, check):
                del capture_output, text, check
                output_dir = temp_path / "out" / "analysis" / "run5"
                output_dir.mkdir(parents=True, exist_ok=True)
                (output_dir / "asv_sequences.fasta").write_text(">ASV1\nACGTACGTACGT\n", encoding="utf-8")
                (output_dir / "read_tracking.csv").write_text("input,filtered,denoisedF,denoisedR,merged,nonchim\n1,1,1,1,1,1\n", encoding="utf-8")
                (output_dir / "asv_table.rds").write_text("stub", encoding="utf-8")
                (output_dir / "asv_abundance.csv").write_text("ASV1\n", encoding="utf-8")
                (output_dir / "taxonomy.csv").write_text("ASV_ID,Taxonomy\nASV1,K;P;C\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, stdout="dada2 ok", stderr="")

            with patch("analysis.dada2_runner.subprocess.run", side_effect=fake_run):
                result = run_dada2([sample], config)

            self.assertTrue(result.success)
            self.assertIn("asv_sequences_fasta", result.files)

    def test_mafft_sanitize_headers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta = temp_path / "input.fasta"
            fasta.write_text(">Sample 1:foo\n" + "A" * 120 + "\n", encoding="utf-8")
            output = temp_path / "sanitized.fasta"
            mapping = sanitize_fasta_for_phylo(fasta, output)
            self.assertTrue(output.exists())
            self.assertIn("ASV1", mapping)
            self.assertEqual(mapping["ASV1"], "Sample 1:foo")

    def test_fasttree_output_exists(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            aligned = temp_path / "aligned.fasta"
            aligned.write_text(">ASV1\n" + "A" * 120 + "\n>ASV2\n" + "T" * 120 + "\n>ASV3\n" + "G" * 120 + "\n>ASV4\n" + "C" * 120 + "\n", encoding="utf-8")

            def fake_run(command, capture_output, text, check):
                del capture_output, text, check
                return subprocess.CompletedProcess(command, 0, stdout="(ASV1,ASV2);", stderr="")

            with patch("phylo.phylo_runner.subprocess.run", side_effect=fake_run):
                result = run_fasttree(aligned, temp_path / "fasttree", {"gtr": True, "gamma": True})

            self.assertTrue(result.success)
            self.assertTrue(result.files["tree_nwk"].exists())
            self.assertGreater(result.files["tree_nwk"].stat().st_size, 0)

    def test_raxml_forbidden_chars(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta = temp_path / "bad.fasta"
            fasta.write_text(">Bad:header,1\nACGT\n", encoding="utf-8")
            errors = validate_fasta_headers(fasta)
            self.assertTrue(errors)

    def test_taxonomy_fields_to_sintax(self):
        result = taxonomy_fields_to_sintax(["Eukaryota", "Chordata", "Actinopteri"])
        self.assertEqual(result, "d:Eukaryota,p:Chordata,c:Actinopteri")

    def test_taxonomy_fields_to_sintax_extends_ranks_beyond_seven(self):
        fields = [f"level{i}" for i in range(1, 9)]
        result = taxonomy_fields_to_sintax(fields)
        self.assertTrue(result.endswith("r8:level8"))
        self.assertEqual(result.count(","), 7)

    def test_bold_skips_if_no_db(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta = temp_path / "asv_sequences.fasta"
            fasta.write_text(">ASV1\n" + "A" * 120 + "\n", encoding="utf-8")
            result = run_boldigger(fasta, {"db": 1, "mode": 2, "thresholds": [97, 95, 90, 85]}, temp_path / "missing.ddb")
            self.assertTrue(result.exists())


if __name__ == "__main__":
    unittest.main()