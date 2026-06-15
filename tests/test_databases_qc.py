"""Contract tests for database and QC runners."""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from core.models import AnalysisMode, Marker, Sample, SamplePaths, StepOutput


class TestFastQC(unittest.TestCase):
	def test_fastqc_se_no_r2_and_failure_is_not_exception(self):
		from qc.fastqc_runner import run_fastqc

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			sample = Sample(
				internal_id="s1",
				display_label="S1",
				paths=SamplePaths(r1=temp_path / "sample_R1.fastq"),
				mode=AnalysisMode.SE,
				marker=Marker.COI,
			)
			with patch("qc.fastqc_runner.subprocess.run") as mocked_run:
				mocked_run.side_effect = FileNotFoundError("fastqc not found")
				result = run_fastqc(sample, temp_path)
			self.assertFalse(result.success)
			self.assertIn("FastQC", result.error_message or "")
			args = mocked_run.call_args[0][0]
			self.assertIn(str(sample.paths.r1), args)
			self.assertNotIn("sample_R2.fastq", " ".join(args))


class TestSILVA(unittest.TestCase):
	def test_silva_missing_file_raises(self):
		from databases.silva_db import load_silva
		with self.assertRaises(ValueError) as ctx:
			load_silva(Path("/nonexistent/silva.fasta"))
		self.assertIn("nao encontrada", str(ctx.exception))


class TestPR2(unittest.TestCase):
	def test_pr2_missing_file_raises(self):
		from databases.pr2_db import load_pr2
		with self.assertRaises(ValueError) as ctx:
			load_pr2(Path("/nonexistent/pr2.fasta"))
		self.assertIn("nao encontrada", str(ctx.exception))


class TestNCBI(unittest.TestCase):
	def test_blast_local_without_db_returns_error(self):
		from databases.ncbi_db import NCBIConfig, run_blast

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			query = temp_path / "query.fasta"
			query.write_text(">ASV1\nACGT\n", encoding="utf-8")
			cfg = NCBIConfig(local_db_path=temp_path / "nt", use_remote=False)
			with patch("databases.ncbi_db.shutil.which", return_value="/usr/bin/blastn"):
				result = run_blast(query, cfg, temp_path)
			self.assertFalse(result.success)
			self.assertIn("update_blastdb", result.error_message or "")


class TestBOLDigger(unittest.TestCase):
	def test_bold_remote_mode_allows_missing_ddb(self):
		from databases.bold_db import BOLDConfig, run_boldigger

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			query = temp_path / "query.fasta"
			query.write_text(">ASV1\nACGT\n", encoding="utf-8")
			bold_db = temp_path / "bold_db"
			bold_db.mkdir()
			cfg = BOLDConfig(db_path=bold_db)

			def fake_run(command, capture_output, text, timeout, cwd, check):
				Path(cwd, "boldigger_result.parquet").write_text("ASV_ID,kingdom\nASV1,Animalia\n", encoding="utf-8")
				Path(cwd, "boldigger_result.xlsx").write_text("stub", encoding="utf-8")
				return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

			with patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"), patch("databases.bold_db.subprocess.run", side_effect=fake_run):
				result = run_boldigger(query, cfg, temp_path)
			self.assertTrue(result.success)
			self.assertIn("bold_parquet", result.files)

	def test_bold_restores_original_headers_after_sanitization(self):
		from databases.bold_db import BOLDConfig, run_boldigger

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			query = temp_path / "query.fasta"
			query.write_text(
				">ASV1;size=13;tax=Animalia\nACGT\n",
				encoding="utf-8",
			)
			bold_db = temp_path / "bold_db"
			bold_db.mkdir()
			cfg = BOLDConfig(db_path=bold_db)

			def fake_run(command, capture_output, text, timeout, cwd, check):
				Path(cwd, "boldigger_result.parquet").write_text("ASV_ID,kingdom\nASV1,Animalia\n", encoding="utf-8")
				return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

			with patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"), patch("databases.bold_db.subprocess.run", side_effect=fake_run):
				result = run_boldigger(query, cfg, temp_path)
			self.assertTrue(result.success)
			remapped = Path(result.files["bold_parquet"]).read_text(encoding="utf-8")
			self.assertIn("ASV1;size=13;tax=Animalia", remapped)


class TestLCAConsensus(unittest.TestCase):
	def test_lca_full_agreement(self):
		from databases.taxonomy_runner import _lca_consensus

		records = [
			{"kingdom": "Animalia", "phylum": "Platyhelminthes", "class_": "Trematoda", "order": "Diplostomida", "family": "Diplostomidae", "genus": "Diplostomum", "species": "Diplostomum spathaceum"},
			{"kingdom": "Animalia", "phylum": "Platyhelminthes", "class_": "Trematoda", "order": "Diplostomida", "family": "Diplostomidae", "genus": "Diplostomum", "species": "Diplostomum spathaceum"},
		]
		result = _lca_consensus(records)
		self.assertEqual(result.get("species"), "Diplostomum spathaceum")

	def test_lca_species_conflict_returns_genus(self):
		from databases.taxonomy_runner import _lca_consensus

		records = [
			{"kingdom": "Animalia", "phylum": "Platyhelminthes", "class_": "Trematoda", "order": "Diplostomida", "family": "Diplostomidae", "genus": "Diplostomum", "species": "Diplostomum spathaceum"},
			{"kingdom": "Animalia", "phylum": "Platyhelminthes", "class_": "Trematoda", "order": "Diplostomida", "family": "Diplostomidae", "genus": "Diplostomum", "species": "Diplostomum pseudospathaceum"},
		]
		result = _lca_consensus(records)
		self.assertNotIn("species", result)
		self.assertEqual(result.get("genus"), "Diplostomum")


class TestTaxonomyPipeline(unittest.TestCase):
	def test_pipeline_allows_zero_sources_and_succeeds_with_mocked_sources(self):
		from databases.taxonomy_runner import run_taxonomy_pipeline

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			fasta = temp_path / "input.fasta"
			fasta.write_text(">ASV1\nACGT\n>ASV2\nTGCA\n", encoding="utf-8")

			silva_csv = temp_path / "silva.csv"
			silva_csv.write_text("ASV_ID,kingdom,genus,species\nASV1,Animalia,Diplostomum,Diplostomum spathaceum\n", encoding="utf-8")
			pr2_csv = temp_path / "pr2.csv"
			pr2_csv.write_text("ASV_ID,kingdom,genus,species\nASV1,Animalia,Diplostomum,Diplostomum spathaceum\n", encoding="utf-8")
			blast_tsv = temp_path / "blast.tsv"
			blast_tsv.write_text("ASV1\thit\t99.0\n", encoding="utf-8")
			bold_parquet = temp_path / "bold.parquet"
			bold_parquet.write_text("ASV_ID,kingdom,genus,species\nASV1,Animalia,Diplostomum,Diplostomum spathaceum\n", encoding="utf-8")
			silva_db = temp_path / "silva.fasta"
			silva_db.write_text(">db\nACGT\n", encoding="utf-8")
			pr2_db = temp_path / "pr2.fasta"
			pr2_db.write_text(">db\nACGT\n", encoding="utf-8")
			ncbi_db = temp_path / "nt"
			(ncbi_db.with_suffix(".nal")).touch()
			bold_db = temp_path / "bold_db"
			bold_db.mkdir()
			(bold_db / "latest.ddb").touch()

			class FakeBold:
				success = True
				files = {"bold_parquet": bold_parquet}
				log = temp_path / "bold.log"
				error_message = None

			class FakeBlast:
				success = True
				files = {"blast_tsv": blast_tsv, "blast_no_hit": temp_path / "nohit.txt", "log": temp_path / "blast.log"}
				log = temp_path / "blast.log"
				error_message = None

			from databases.bold_db import BOLDConfig
			from databases.ncbi_db import NCBIConfig
			from databases.pr2_db import PR2Config
			from databases.silva_db import SILVAConfig

			silva_cfg = SILVAConfig(db_path=silva_db)
			pr2_cfg = PR2Config(db_path=pr2_db)
			ncbi_cfg = NCBIConfig(local_db_path=ncbi_db, use_remote=False)
			bold_cfg = BOLDConfig(db_path=bold_db)

			with patch("databases.taxonomy_runner.run_sintax", side_effect=[silva_csv, pr2_csv]), \
				 patch("databases.taxonomy_runner.run_blast", return_value=FakeBlast()), \
				 patch("databases.taxonomy_runner.run_boldigger_step", return_value=FakeBold()):
				result = run_taxonomy_pipeline(fasta, temp_path / "taxonomy", silva_config=silva_cfg, pr2_config=pr2_cfg, ncbi_config=ncbi_cfg, bold_config=bold_cfg)

			self.assertTrue(result.success)
			self.assertIn("taxonomy_consensus_csv", result.files)
			self.assertIn("annotated_fasta", result.files)
			self.assertIn("conflict_report_csv", result.files)

	def test_pipeline_without_sources_fails_cleanly(self):
		from databases.taxonomy_runner import run_taxonomy_pipeline

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			fasta = temp_path / "input.fasta"
			fasta.write_text(">ASV1\nACGT\n", encoding="utf-8")
			result = run_taxonomy_pipeline(fasta, temp_path / "taxonomy")
			self.assertFalse(result.success)
			self.assertIn("Nenhuma fonte", result.error_message or "")


if __name__ == "__main__":
	unittest.main()