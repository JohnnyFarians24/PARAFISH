"""Phase 3 regression tests: download path traversal, BOLD credentials, phylo keys."""

from __future__ import annotations

import subprocess
import tempfile
import unittest

import pytest
from pathlib import Path
from unittest.mock import patch

from core.models import StepOutput


class DownloadTraversalTests(unittest.TestCase):
	def test_run_zip_rejects_paths_outside_runs_root(self):
		with tempfile.TemporaryDirectory() as temp_dir:
			runs = Path(temp_dir) / "runs"
			runs.mkdir()
			with patch.dict("os.environ", {"RUNS_DIR": str(runs)}):
				from app import create_app

				client = create_app().server.test_client()
				self.assertEqual(client.get("/download/run_zip?run_dir=/etc").status_code, 403)
				self.assertEqual(client.get(f"/download/run_zip?run_dir={runs}/../secrets").status_code, 403)
				self.assertEqual(client.get("/download/run_zip").status_code, 400)

	def test_runs_route_rejects_traversal(self):
		with tempfile.TemporaryDirectory() as temp_dir:
			runs = Path(temp_dir) / "runs"
			runs.mkdir()
			with patch.dict("os.environ", {"RUNS_DIR": str(runs)}):
				from app import create_app

				client = create_app().server.test_client()
				self.assertEqual(client.get("/runs/../../etc/passwd").status_code, 403)

	def test_run_zip_allows_legit_run_dir(self):
		with tempfile.TemporaryDirectory() as temp_dir:
			runs = Path(temp_dir) / "runs"
			run = runs / "run_ok"
			run.mkdir(parents=True)
			(run / "progress.json").write_text("{}", encoding="utf-8")
			with patch.dict("os.environ", {"RUNS_DIR": str(runs)}):
				from app import create_app

				client = create_app().server.test_client()
				response = client.get(f"/download/run_zip?run_dir={run}")
				self.assertEqual(response.status_code, 200)
				self.assertEqual(response.mimetype, "application/zip")


class BOLDCredentialTests(unittest.TestCase):
	def test_bold_identify_runs_without_credentials(self):
		"""``identify`` queries BOLD without logging in; only ``download_db`` needs an account."""
		from databases.bold_db import BOLDConfig, run_boldigger

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			query = temp_path / "query.fasta"
			query.write_text(">ASV1\nACGT\n", encoding="utf-8")
			bold_db = temp_path / "bold_db"
			bold_db.mkdir()
			out_dir = temp_path / "out"

			def fake_run(command, **kwargs):
				cwd = kwargs.get("cwd")
				if cwd:
					Path(cwd, "sanitized_input_identification_result.parquet.snappy").write_text(
						"ASV_ID,kingdom\nASV1,Animalia\n", encoding="utf-8"
					)
				return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

			with patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"), \
				patch("databases.bold_db.boldigger_version", return_value=None), \
				patch("databases.bold_db.get_bold_credentials", return_value=None), \
				patch("databases.bold_db.subprocess.run", side_effect=fake_run) as mocked_run:
				result = run_boldigger(query, BOLDConfig(db_path=bold_db), out_dir)

			self.assertTrue(result.success, result.error_message)
			mocked_run.assert_called_once()
			self.assertIn("bold_parquet", result.files)

	def test_bold_rejects_versions_with_upstream_crash_bug(self):
		from databases.bold_db import BOLDConfig

		with tempfile.TemporaryDirectory() as temp_dir:
			with patch("databases.bold_db.boldigger_available", return_value=True), \
				patch("databases.bold_db.boldigger_version", return_value=(1, 6, 0)):
				errors = BOLDConfig(db_path=Path(temp_dir)).validate()
		self.assertTrue(any("pbar.total" in error for error in errors), errors)

	def test_bold_does_not_require_a_populated_local_ddb(self):
		"""A missing/empty local DB must not block BOLD.

		In boldigger3 v3 the local DuckDB feeds ONLY ``add_metadata`` (sex,
		institution, country). Identification — taxonomy and % identity, the only
		things this pipeline consumes — comes from BOLD's remote API. Refusing to run
		without the 2.3 GB database therefore blocked runs that work perfectly well,
		and on a small machine that database cannot even be built (the OOM killer
		takes it out mid-build, leaving an empty .ddb behind).
		"""
		from databases.bold_db import BOLDConfig

		with tempfile.TemporaryDirectory() as temp_dir:
			with patch("databases.bold_db.boldigger_available", return_value=True), \
				patch("databases.bold_db.boldigger_version", return_value=(3, 0, 2)):
				errors = BOLDConfig(db_path=Path(temp_dir)).validate()
		self.assertEqual(errors, [])

	def test_empty_ddb_is_repaired_with_a_schema_only_metadata_table(self):
		"""An empty .ddb gets a well-formed, zero-row ``bold_public``.

		Without the table, ``add_metadata``'s LEFT JOIN aborts with "table bold_public
		does not exist" and destroys an identification that already succeeded. With an
		empty one, the join yields NULL metadata and the run completes.
		"""
		duckdb = pytest.importorskip("duckdb")
		from databases.bold_db import bold_metadata_db_status, ensure_bold_metadata_db

		with tempfile.TemporaryDirectory() as temp_dir:
			db_dir = Path(temp_dir)
			empty = db_dir / "BOLD_Public.ddb"
			duckdb.connect(str(empty)).close()  # a table-less DuckDB, as the OOM leaves it
			self.assertEqual(bold_metadata_db_status(empty), (False, 0))

			path, note = ensure_bold_metadata_db(db_dir)
			self.assertIsNotNone(path)
			self.assertIn("is empty", (note or ""))
			self.assertEqual(bold_metadata_db_status(path), (True, 0))

	def test_boldigger_command_passes_db_path_on_modern_boldigger(self):
		from databases import bold_db
		from databases.bold_db import BOLDConfig, _boldigger_command

		with tempfile.TemporaryDirectory() as temp_dir:
			ddb = Path(temp_dir) / "bold.ddb"
			ddb.write_text("", encoding="utf-8")
			# Pin BOLDIGGER_PYTHON off: it is read from the developer's .env and now
			# takes precedence over PATH, which would otherwise leak into this assertion.
			with patch.object(bold_db, "BOLDIGGER_PYTHON", ""), \
				patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"), \
				patch("databases.bold_db.boldigger_version", return_value=(3, 0, 2)):
				command = _boldigger_command(Path("/tmp/in.fasta"), BOLDConfig(db_path=Path(temp_dir)))
		# identify <fasta> <db_path> --db N --mode M --thresholds ...
		self.assertEqual(command[:4], ["boldigger3", "identify", "/tmp/in.fasta", str(ddb)])
		self.assertIn("--db", command)

	def test_boldigger_command_has_no_extra_positional_on_legacy(self):
		from databases import bold_db
		from databases.bold_db import BOLDConfig, _boldigger_command

		with patch.object(bold_db, "BOLDIGGER_PYTHON", ""), \
			patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"), \
			patch("databases.bold_db.boldigger_version", return_value=(1, 6, 0)):
			command = _boldigger_command(Path("/tmp/in.fasta"), BOLDConfig(db_path=Path("/tmp/db")))
		# The 1.x CLI takes exactly one positional (the FASTA).
		self.assertEqual(command[:3], ["boldigger3", "identify", "/tmp/in.fasta"])
		self.assertNotIn("/tmp/db", command)
		self.assertIn("--db", command)


class PhyloKeyNormalizationTests(unittest.TestCase):
	def test_raxml_aliases_tree_nwk_to_support_tree(self):
		from phylo.phylo_runner import _finalize_phylo_output

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			support = temp_path / "raxml.raxml.supportFBP"
			support.write_text("(A,B);", encoding="utf-8")
			best = temp_path / "raxml.raxml.bestTree"
			best.write_text("(A,B);", encoding="utf-8")
			out = _finalize_phylo_output(
				StepOutput(step_name="raxml", success=True, files={"best_tree": best, "support_tree": support})
			)
			self.assertEqual(out.files["tree_nwk"], support)

	def test_raxml_falls_back_to_best_tree_when_support_missing(self):
		from phylo.phylo_runner import _finalize_phylo_output

		with tempfile.TemporaryDirectory() as temp_dir:
			temp_path = Path(temp_dir)
			best = temp_path / "raxml.raxml.bestTree"
			best.write_text("(A,B);", encoding="utf-8")
			out = _finalize_phylo_output(
				StepOutput(step_name="raxml", success=True, files={"best_tree": best, "support_tree": temp_path / "missing"})
			)
			self.assertEqual(out.files["tree_nwk"], best)

	def test_fasttree_tree_nwk_preserved(self):
		from phylo.phylo_runner import _finalize_phylo_output

		with tempfile.TemporaryDirectory() as temp_dir:
			nwk = Path(temp_dir) / "tree.nwk"
			out = _finalize_phylo_output(
				StepOutput(step_name="fasttree", success=True, files={"tree_nwk": nwk})
			)
			self.assertEqual(out.files["tree_nwk"], nwk)



class BoldiggerVersionDetectionTests(unittest.TestCase):
    def _version_from(self, stdout: str, stderr: str = ""):
        from databases import bold_db

        completed = subprocess.CompletedProcess(["boldigger3", "--version"], 0, stdout=stdout, stderr=stderr)
        with patch("databases.bold_db.boldigger_available", return_value=True), \
            patch("databases.bold_db.subprocess.run", return_value=completed):
            bold_db.boldigger_version.cache_clear()
            try:
                return bold_db.boldigger_version()
            finally:
                bold_db.boldigger_version.cache_clear()

    def test_version_is_read_from_its_own_line(self):
        noisy = "Installing dependencies...\nYour boldigger3 version is outdated.\n1.6.0\n"
        self.assertEqual(self._version_from(noisy), (1, 6, 0))

    def test_version_ignores_python_version_inside_traceback_paths(self):
        """A loose search picks "3.12" out of .../lib/python3.12/site-packages/..."""
        traceback_text = (
            'Traceback (most recent call last):\n'
            '  File "/opt/venv/lib/python3.12/site-packages/boldigger3/__main__.py", line 2\n'
            'ImportError: lxml.html.clean module is now a separate project\n'
        )
        self.assertIsNone(self._version_from("", traceback_text))

    def test_modern_version_is_detected(self):
        self.assertEqual(self._version_from("3.0.2\n"), (3, 0, 2))


class BoldiggerInterpreterPrecedenceTests(unittest.TestCase):
    def test_explicit_interpreter_wins_over_path(self):
        """BOLDIGGER_PYTHON exists to point at another env; a stale boldigger3 on
        PATH must not shadow it."""
        from databases import bold_db

        with patch.object(bold_db, "BOLDIGGER_PYTHON", "/opt/boldenv/bin/python"), \
            patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"):
            self.assertEqual(bold_db._boldigger_base_command(), ["/opt/boldenv/bin/python", "-m", "boldigger3"])

    def test_path_is_used_when_no_interpreter_configured(self):
        from databases import bold_db

        with patch.object(bold_db, "BOLDIGGER_PYTHON", ""), \
            patch("databases.bold_db.shutil.which", return_value="/usr/bin/boldigger3"):
            self.assertEqual(bold_db._boldigger_base_command(), ["boldigger3"])

if __name__ == "__main__":
	unittest.main()
