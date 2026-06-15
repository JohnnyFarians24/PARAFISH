"""Contract tests for the phase 2 core data model and validation helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from core.models import AnalysisBranch, AnalysisMode, Marker, PhyloMode, RunConfig, Sample, SamplePaths
from core.validation import validate_fasta_headers, validate_manifest, validate_sample_names


class Phase2ContractTests(unittest.TestCase):
    def _make_fastq(self, directory: Path, name: str, reads: int = 1) -> Path:
        path = directory / name
        lines = []
        for index in range(reads):
            lines.extend([
                f"@read{index + 1}",
                "ACGT",
                "+",
                "!!!!",
            ])
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def test_runconfig_vsearch_pe(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            r1 = self._make_fastq(temp_path, "sample_R1.fastq", reads=2)
            r2 = self._make_fastq(temp_path, "sample_R2.fastq", reads=2)
            sample = Sample(
                internal_id="sample_001",
                display_label="Sample 1",
                paths=SamplePaths(r1=r1, r2=r2),
                mode=AnalysisMode.PE,
                marker=Marker.COI,
            )
            config = RunConfig(
                run_id="run_001",
                samples=[sample],
                branch=AnalysisBranch.VSEARCH,
                phylo_mode=PhyloMode.FAST,
                output_dir=temp_path / "out",
            )

            self.assertEqual(config.branch, AnalysisBranch.VSEARCH)
            self.assertEqual(config.samples[0].mode, AnalysisMode.PE)
            self.assertEqual(config.samples[0].marker, Marker.COI)

    def test_runconfig_dada2_pe(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            r1 = self._make_fastq(temp_path, "sample_R1.fastq", reads=2)
            r2 = self._make_fastq(temp_path, "sample_R2.fastq", reads=2)
            ref_fasta = temp_path / "ref.fasta"
            ref_fasta.write_text(">ref\nACGT\n", encoding="utf-8")
            sample = Sample(
                internal_id="sample_002",
                display_label="Sample 2",
                paths=SamplePaths(r1=r1, r2=r2),
                mode=AnalysisMode.PE,
                marker=Marker.SSU_18S,
            )
            config = RunConfig(
                run_id="run_002",
                samples=[sample],
                branch=AnalysisBranch.DADA2,
                phylo_mode=PhyloMode.ROBUST,
                output_dir=temp_path / "out",
            )
            config.dada2_params["ref_fasta"] = ref_fasta

            self.assertEqual(config.branch, AnalysisBranch.DADA2)
            self.assertEqual(config.phylo_mode, PhyloMode.ROBUST)
            self.assertEqual(config.dada2_params["ref_fasta"], ref_fasta)

    def test_validate_manifest_missing_r2_in_pe(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            r1 = self._make_fastq(temp_path, "sample_R1.fastq", reads=2)
            sample = Sample(
                internal_id="sample_003",
                display_label="Sample 3",
                paths=SamplePaths(r1=r1, r2=None),
                mode=AnalysisMode.PE,
                marker=Marker.COI,
            )
            config = RunConfig(
                run_id="run_003",
                samples=[sample],
                branch=AnalysisBranch.VSEARCH,
                output_dir=temp_path / "out",
            )

            errors = validate_manifest(config)
            self.assertTrue(any("Modo PE requer R2" in error or "R2" in error for error in errors))

    def test_validate_fasta_headers_detects_colon(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta = temp_path / "test.fasta"
            fasta.write_text(">Sample:1\nACGT\n>Sample 2\nTGCA\n", encoding="utf-8")

            errors = validate_fasta_headers(fasta)
            self.assertIn("Sample:1", errors)

    def test_validate_sample_names_rejects_bad_internal_ids(self):
        samples = [
            Sample(
                internal_id="bad id",
                display_label="Bad",
                paths=SamplePaths(r1=Path("/tmp/missing.fastq")),
                mode=AnalysisMode.SE,
                marker=Marker.OTHER,
            )
        ]
        errors = validate_sample_names(samples)
        self.assertTrue(errors)


if __name__ == "__main__":
    unittest.main()