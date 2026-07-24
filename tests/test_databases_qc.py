"""Contract tests for database and QC runners."""

from __future__ import annotations

import csv
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

	def _refusal_run(self, stderr: str, temp_path: Path):
		"""run_blast against a blastn that EXITS 0 while printing ``stderr``."""
		import subprocess

		from databases.ncbi_db import NCBIConfig, run_blast

		query = temp_path / "query.fasta"
		query.write_text(">OTU1\nACGT\n", encoding="utf-8")
		completed = subprocess.CompletedProcess(args=["blastn"], returncode=0, stdout="", stderr=stderr)
		cfg = NCBIConfig(use_remote=True)
		with patch("databases.ncbi_db.shutil.which", return_value="/usr/bin/blastn"), \
			patch("databases.ncbi_db.subprocess.run", return_value=completed), \
			patch("databases.ncbi_db.time.sleep"), \
			patch("databases.ncbi_db.REMOTE_BACKOFF_SEC", (0, 0)):
			return run_blast(query, cfg, temp_path)

	def test_a_server_refusal_is_a_failure_even_though_blastn_exits_zero(self):
		"""blastn exits 0 having printed "Severe Error ... SIGXCPU" for every query and
		written an empty result. Read only the exit code and the step reports success with
		hit_count=0 — and NCBI's refusal reaches source_status wearing the word "no_hits",
		which is the exact lie this layer exists to prevent."""
		refusal = (
			"Error: [blastn] internal_error: (Severe Error) Blast search error: Details: search failed. "
			"# Informational Message: [blastsrv4.REAL]: Error: CPU usage limit was exceeded, "
			"resulting in SIGXCPU (24)."
		)
		with tempfile.TemporaryDirectory() as temp_dir:
			result = self._refusal_run(refusal, Path(temp_dir))
		self.assertFalse(result.success)
		self.assertIn("recusou", (result.error_message or "").lower())
		self.assertIn("megablast", (result.error_message or "").lower())

	def test_a_clean_empty_result_is_not_treated_as_a_refusal(self):
		"""The mirror case: blastn ran fine and this query genuinely matches nothing in nt.
		That IS an absence of hits, and must stay a success with hit_count=0."""
		with tempfile.TemporaryDirectory() as temp_dir:
			result = self._refusal_run("", Path(temp_dir))
		self.assertTrue(result.success)
		self.assertEqual(result.metrics.get("hit_count"), 0.0)

	def test_remote_default_asks_for_megablast_and_no_word_size(self):
		"""word_size 11 against all of nt is a full blastn search — the most expensive
		thing you can ask a shared public quota for, and it is refused (SIGXCPU) while
		megablast on the same query returns 123 hits."""
		from databases.ncbi_db import NCBIConfig

		cfg = NCBIConfig(use_remote=True)
		self.assertEqual(cfg.task, "megablast")
		self.assertIsNone(cfg.word_size)


class TestBOLDigger(unittest.TestCase):
	def test_bold_remote_mode_allows_missing_ddb(self):
		"""A legacy (remote-API) boldigger3 must not be blocked on a local *.ddb."""
		from databases.bold_db import BOLDConfig, run_boldigger
		from core.config import BOLDCredentials

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			query = temp_path / "query.fasta"
			query.write_text(">ASV1\nACGT\n", encoding="utf-8")
			bold_db = temp_path / "bold_db"
			bold_db.mkdir()
			cfg = BOLDConfig(db_path=bold_db)

			def fake_run(command, **kwargs):
				cwd = kwargs["cwd"]
				Path(cwd, "sanitized_input_identification_result.parquet.snappy").write_text("ASV_ID,kingdom\nASV1,Animalia\n", encoding="utf-8")
				Path(cwd, "sanitized_input_identification_result.xlsx").write_text("stub", encoding="utf-8")
				return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

			with patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"), \
				patch("databases.bold_db.boldigger_version", return_value=None), \
				patch("databases.bold_db.get_bold_credentials", return_value=BOLDCredentials(username="u", password="p")), \
				patch("databases.bold_db.subprocess.run", side_effect=fake_run):
				result = run_boldigger(query, cfg, temp_path)
			self.assertTrue(result.success, result.error_message)
			self.assertIn("bold_parquet", result.files)

	def test_bold_finds_results_in_boldigger3_data_subdir(self):
		"""boldigger3 3.x writes its results into a ``boldigger3_data`` subdirectory."""
		from databases.bold_db import BOLDConfig, run_boldigger

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			query = temp_path / "query.fasta"
			query.write_text(">ASV1\nACGT\n", encoding="utf-8")
			bold_db = temp_path / "bold_db"
			bold_db.mkdir()
			(bold_db / "public.ddb").write_text("", encoding="utf-8")

			def fake_run(command, **kwargs):
				data_dir = Path(kwargs["cwd"]) / "boldigger3_data"
				data_dir.mkdir(parents=True, exist_ok=True)
				(data_dir / "sanitized_input_identification_result.parquet.snappy").write_text("ASV_ID,kingdom\nASV1,Animalia\n", encoding="utf-8")
				return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

			with patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"), \
				patch("databases.bold_db.boldigger_version", return_value=(3, 0, 2)), \
				patch("databases.bold_db.get_bold_credentials", return_value=None), \
				patch("databases.bold_db.subprocess.run", side_effect=fake_run):
				result = run_boldigger(query, BOLDConfig(db_path=bold_db), temp_path / "out")
			self.assertTrue(result.success, result.error_message)
			self.assertEqual(Path(result.files["bold_parquet"]).parent.name, "boldigger3_data")

	def test_bold_restores_original_headers_after_sanitization(self):
		from databases.bold_db import BOLDConfig, run_boldigger
		from core.config import BOLDCredentials

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

			def fake_run(command, **kwargs):
				Path(kwargs["cwd"], "sanitized_input_identification_result.parquet.snappy").write_text("ASV_ID,kingdom\nASV1,Animalia\n", encoding="utf-8")
				return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

			with patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"), \
				patch("databases.bold_db.boldigger_version", return_value=None), \
				patch("databases.bold_db.get_bold_credentials", return_value=BOLDCredentials(username="u", password="p")), \
				patch("databases.bold_db.subprocess.run", side_effect=fake_run):
				result = run_boldigger(query, cfg, temp_path)
			self.assertTrue(result.success, result.error_message)
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

class SourceStatusTruthTests(unittest.TestCase):
    """A source may not report success it did not deliver.

    Every source used to judge itself by a different proxy: BOLD by whether its result
    FILE existed, NCBI by BLAST's hit_count, custom BLAST by a match counter. So BOLD
    reported "ok, identificação BOLD concluída" on a run where all 7 of its rows were
    unassigned with identity 0.0.
    """

    def _records(self, **by_source):
        # {source: [asv_ids...]} -> records_by_asv, with a real taxon on each record
        records = {}
        for source, asv_ids in by_source.items():
            for asv_id in asv_ids:
                records.setdefault(asv_id, []).append({"source": source, "genus": "Gadus", "species": "Gadus morhua"})
        return records

    def test_a_source_with_zero_assignments_cannot_stay_ok(self):
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "bold.parquet"
            raw.write_text("x" * 100, encoding="utf-8")  # BOLD did return a file...
            status = [{"source": "BOLD", "status": "ok", "detail": "identificação BOLD concluída", "key": "BOLD", "raw": raw}]
            out = _reconcile_source_status(status, self._records(SILVA=["OTU1"]), 7)

        self.assertEqual(out[0]["status"], "no_hits")  # ...but identified nothing
        self.assertIn("7", out[0]["detail"])

    def test_a_source_arriving_as_no_hits_with_empty_output_is_downgraded_to_failed(self):
        """Reconciliation used to judge only sources claiming "ok". NCBI labels itself
        "ok if hits > 0 else no_hits" up front, so a remote BLAST killed server-side
        ("CPU usage limit was exceeded" — 0-byte result, error-filled log) arrived here
        already saying no_hits and was waved through: a failure wearing the words of an
        honest absence."""
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "blast_result.tsv"
            raw.write_text("", encoding="utf-8")  # BLAST died; nothing came back
            status = [{"source": "NCBI", "status": "no_hits", "detail": "0 hits BLAST", "key": "NCBI", "raw": raw}]
            out = _reconcile_source_status(status, self._records(SILVA=["OTU1"]), 7, marker="18S")

        self.assertEqual(out[0]["status"], "failed")

    def test_a_source_arriving_as_no_hits_with_real_output_stays_no_hits(self):
        """The mirror case: real output and no assignments IS an honest absence."""
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "blast_result.tsv"
            raw.write_text("some real blast output", encoding="utf-8")
            status = [{"source": "NCBI", "status": "no_hits", "detail": "0 hits", "key": "NCBI", "raw": raw}]
            out = _reconcile_source_status(status, self._records(SILVA=["OTU1"]), 7, marker="18S")

        self.assertEqual(out[0]["status"], "no_hits")

    def test_bold_on_18s_is_not_applicable_not_no_hits(self):
        """BOLD indexes COI. On 18S its 0/7 is not an absence of matches — it is a
        question it was never able to answer, and saying "no_hits" invites the reader
        to conclude the taxa are absent from BOLD."""
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "bold.parquet"
            raw.write_text("x" * 100, encoding="utf-8")
            status = [{"source": "BOLD", "status": "ok", "detail": "concluída", "key": "BOLD", "raw": raw}]
            out = _reconcile_source_status(status, self._records(SILVA=["OTU1"]), 7, marker="18S")

        self.assertEqual(out[0]["status"], "not_applicable")
        self.assertIn("COI", out[0]["detail"])
        self.assertIn("18S", out[0]["detail"])

    def test_bold_on_coi_with_no_assignments_is_still_no_hits(self):
        """The marker matches, so 0/7 IS an honest absence and must not be excused."""
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "bold.parquet"
            raw.write_text("x" * 100, encoding="utf-8")
            status = [{"source": "BOLD", "status": "ok", "detail": "concluída", "key": "BOLD", "raw": raw}]
            out = _reconcile_source_status(status, self._records(SILVA=["OTU1"]), 7, marker="COI")

        self.assertEqual(out[0]["status"], "no_hits")

    def test_an_empty_raw_output_is_a_failure_not_an_absence_of_hits(self):
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "custom_blast_raw.tsv"
            raw.write_text("", encoding="utf-8")  # 0 bytes: BLAST produced nothing
            status = [{"source": "Base customizada (BLAST)", "status": "ok", "detail": "0/7", "key": "CUSTOM", "raw": raw}]
            out = _reconcile_source_status(status, {}, 7)

        self.assertEqual(out[0]["status"], "failed")
        self.assertIn("vazia", out[0]["detail"])

    def test_a_missing_raw_output_is_also_a_failure(self):
        from databases.taxonomy_runner import _reconcile_source_status

        status = [{"source": "NCBI", "status": "ok", "detail": "68 hits BLAST", "key": "NCBI", "raw": Path("/nonexistent.tsv")}]
        out = _reconcile_source_status(status, {}, 7)

        self.assertEqual(out[0]["status"], "failed")

    def test_a_source_that_really_assigned_keeps_ok_and_gets_a_count(self):
        from databases.taxonomy_runner import _reconcile_source_status

        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "pr2.csv"
            raw.write_text("x", encoding="utf-8")
            status = [{"source": "PR2", "status": "ok", "detail": "7 sequências classificadas", "key": "PR2", "raw": raw}]
            out = _reconcile_source_status(status, self._records(PR2=["OTU1", "OTU2"]), 7)

        self.assertEqual(out[0]["status"], "ok")
        self.assertIn("2/7", out[0]["detail"])

    def test_rows_with_every_rank_empty_are_not_assignments(self):
        from databases.taxonomy_runner import _source_assignment_counts

        records = {"OTU1": [{"source": "BOLD", "genus": "", "species": ""}, {"source": "PR2", "species": "Gadus morhua"}]}
        counts = _source_assignment_counts(records)

        self.assertNotIn("BOLD", counts)
        self.assertEqual(counts["PR2"], 1)

    def test_a_skipped_source_is_left_alone(self):
        # SILVA skipped for RAM must keep saying so, not be relabelled "no_hits".
        from databases.taxonomy_runner import _reconcile_source_status

        status = [{"source": "SILVA", "status": "skipped", "detail": "ignorada: base maior que a RAM"}]
        out = _reconcile_source_status(status, {}, 7)

        self.assertEqual(out[0]["status"], "skipped")
        self.assertIn("ignorada", out[0]["detail"])

    def test_reconciled_rows_carry_only_the_csv_columns(self):
        # key/raw are internal; DictWriter raises on extras.
        from databases.taxonomy_runner import _reconcile_source_status

        status = [{"source": "PR2", "status": "ok", "detail": "x", "key": "PR2", "raw": None}]
        out = _reconcile_source_status(status, {}, 7)

        self.assertEqual(set(out[0]), {"source", "status", "detail"})


class MarkerResolutionTests(unittest.TestCase):
    """Task 1(a): the run's marker comes from the primer preset, not the file name.

    The COI/Leray run that motivated this had samples named ``aptl64b_01``. No
    ``coi`` token in the stem, so the name-based guess returned None, the marker
    resolved to OTHER, and SILVA/PR2 were cut with COI primers.
    """

    def test_preset_supplies_marker_when_file_names_are_silent(self):
        from databases.reference_filter import resolve_run_marker

        from analysis.sample_parsing import parse_fastq_name

        self.assertIsNone(parse_fastq_name("aptl64b_01.fastq.gz").marker)
        marker, note = resolve_run_marker([None], {"primer_preset": "COI_Leray"})

        self.assertEqual(marker, "COI")
        self.assertEqual(note, "")

    def test_preset_wins_over_file_name_and_records_the_disagreement(self):
        from databases.reference_filter import resolve_run_marker

        marker, note = resolve_run_marker(["18S"], {"primer_preset": "COI_Leray"})

        self.assertEqual(marker, "COI")
        self.assertIn("COI", note)
        self.assertIn("18S", note)

    def test_file_name_still_used_when_no_preset(self):
        from databases.reference_filter import resolve_run_marker

        self.assertEqual(resolve_run_marker(["18S"], {})[0], "18S")

    def test_every_preset_declares_a_marker(self):
        from analysis.vsearch_runner import PRESET_MARKER, PRIMER_PRESETS

        self.assertEqual(set(PRESET_MARKER), set(PRIMER_PRESETS))

    def test_unknown_marker_does_not_cut_a_reference_with_arbitrary_primers(self):
        """OTHER + COI primers used to reach the amplicon cut and shred SILVA."""
        from analysis.vsearch_runner import PRIMER_PRESETS
        from databases.reference_filter import build_marker_reference, marker_filter_applies

        leray = PRIMER_PRESETS["COI_Leray"]
        self.assertFalse(marker_filter_applies("OTHER", leray))
        self.assertFalse(marker_filter_applies("", leray))
        with tempfile.TemporaryDirectory() as d:
            source = Path(d) / "silva.fasta"
            source.write_text(">a;tax=d:Eukaryota;\nACGT\n", encoding="utf-8")
            self.assertIsNone(build_marker_reference(source, "OTHER", leray))


class MarkerApplicabilityTests(unittest.TestCase):
    """Task 1(b): applicability is symmetric, and decided before any work."""

    def test_rrna_sources_are_not_applicable_on_a_coi_run(self):
        from databases.taxonomy_runner import _marker_incompatible

        self.assertTrue(_marker_incompatible("SILVA", "COI"))
        self.assertTrue(_marker_incompatible("PR2", "COI"))

    def test_coi_source_is_not_applicable_on_an_rrna_run(self):
        """The direction that already worked must keep working."""
        from databases.taxonomy_runner import _marker_incompatible

        self.assertTrue(_marker_incompatible("BOLD", "18S"))

    def test_each_source_is_applicable_to_its_own_marker(self):
        from databases.taxonomy_runner import _marker_incompatible

        self.assertFalse(_marker_incompatible("SILVA", "18S"))
        self.assertFalse(_marker_incompatible("SILVA", "16S"))
        self.assertFalse(_marker_incompatible("PR2", "18S"))
        self.assertFalse(_marker_incompatible("BOLD", "COI"))

    def test_unknown_marker_never_makes_a_source_inapplicable(self):
        """A gap in configuration must not silently drop a database."""
        from databases.taxonomy_runner import _marker_incompatible

        for marker in ("", "OTHER", None):
            self.assertFalse(_marker_incompatible("SILVA", marker))
            self.assertFalse(_marker_incompatible("BOLD", marker))

    def test_incompatible_source_is_skipped_without_being_filtered_or_loaded(self):
        """The whole point: not_applicable costs nothing. If sintax_ready_path or
        _apply_marker_filter is reached on a COI run, the minutes are already spent."""
        from unittest.mock import MagicMock

        import databases.taxonomy_runner as tr

        silva = MagicMock()
        with tempfile.TemporaryDirectory() as d:
            fasta = Path(d) / "reps.fasta"
            fasta.write_text(">OTU1\nACGTACGTACGT\n", encoding="utf-8")
            with patch.object(tr, "_apply_marker_filter") as filt, patch.object(tr, "run_sintax") as sintax:
                result = tr.run_taxonomy_pipeline(
                    fasta, Path(d) / "out", silva_config=silva, marker="COI"
                )

                status_csv = result.files["source_status_csv"]
                rows = list(csv.DictReader(status_csv.open(encoding="utf-8")))

        silva.sintax_ready_path.assert_not_called()
        filt.assert_not_called()
        sintax.assert_not_called()
        entry = [row for row in rows if row["source"] == "SILVA"][0]
        self.assertEqual(entry["status"], "not_applicable")
        self.assertIn("COI", entry["detail"])

    def test_preflight_announces_the_skip_instead_of_promising_a_filter(self):
        """The reported symptom was a preflight line promising to reduce SILVA "ao teu
        marcador (OTHER)". On a COI run SILVA is not filtered at all — it is skipped."""
        from core.preflight import preflight_notes

        with patch("core.preflight.available_memory_bytes", return_value=2 * 1024**3), \
            patch("core.preflight.total_memory_bytes", return_value=3 * 1024**3):
            notes = preflight_notes(
                ["silva", "pr2", "bold"], object(), 1, "vsearch",
                marker_filter="COI", inapplicable=["SILVA", "PR2"],
            )

        joined = " ".join(notes)
        self.assertIn("SILVA", joined)
        self.assertIn("não aplicável", joined)
        self.assertNotIn("reduzida ao teu marcador", joined)


class SetupMarkerReadoutTests(unittest.TestCase):
    """Addendum to Task 1: the Setup page's marker readouts agree with the run.

    There is no marker *control* on Setup — only displays. They showed the file-name
    guess, so a COI/Leray run read "NONE" on screen while the pipeline acted on OTHER.
    """

    def test_marker_column_never_renders_the_literal_none(self):
        from ui.components.sample_table import _marker_badge

        self.assertEqual(_marker_badge(None).children, "OTHER")
        self.assertEqual(_marker_badge("").children, "OTHER")

    def test_marker_column_shows_the_resolved_run_marker(self):
        from ui.components.sample_table import render_sample_table

        samples = [{"internal_id": "aptl64b_01", "marker": None, "selected": True}]
        rendered = str(render_sample_table(samples, selectable=True, run_marker="COI"))

        self.assertIn("COI", rendered)
        self.assertNotIn("NONE", rendered)

    def test_marker_column_falls_back_to_the_file_name_guess(self):
        """The Samples page has no primer preset to resolve from."""
        from ui.components.sample_table import render_sample_table

        samples = [{"internal_id": "s_18S_01", "marker": "18S", "selected": True}]
        self.assertIn("18S", str(render_sample_table(samples, selectable=True)))
