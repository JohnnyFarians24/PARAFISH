"""BOLD reference database helpers."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

try:
	import pandas as pd
except Exception:  # pragma: no cover
	pd = None

from core.config import BOLD_CACHE_ROOT, BOLDIGGER_PYTHON, BOLDIGGER_TIMEOUT_SEC
from core.models import StepOutput


@dataclass(frozen=True)
class BOLDConfig:
	db_path: Path
	db_source: int = 1
	mode: int = 2
	thresholds: List[int] = field(default_factory=lambda: [97, 95, 90, 85])
	timeout: int = BOLDIGGER_TIMEOUT_SEC

	def validate(self) -> List[str]:
		errors: List[str] = []
		boldigger_available = bool(shutil.which("boldigger3"))
		python_available = bool(BOLDIGGER_PYTHON and (shutil.which(BOLDIGGER_PYTHON) or os.path.isfile(BOLDIGGER_PYTHON)))
		if not boldigger_available and not python_available:
			errors.append("boldigger3 nao encontrado no PATH e BOLDIGGER_PYTHON nao esta definido. Instala com: pip install boldigger3")
		if self.db_source not in range(1, 9):
			errors.append(f"db_source deve ser 1-8, recebido: {self.db_source}")
		if self.mode not in (1, 2, 3):
			errors.append(f"mode deve ser 1, 2 ou 3, recebido: {self.mode}")
		if len(self.thresholds) != 4:
			errors.append(f"thresholds deve ter 4 valores, recebido: {self.thresholds}")
		return errors


def _resolve_db_argument(config: BOLDConfig) -> Path:
	# Always pass the configured db_path directly to boldigger3. The
	# installed boldigger3 (3.x) uses remote API and does not require local
	# DuckDB files; the CLI still expects a positional path argument.
	return config.db_path


def _sanitize_fasta_for_bold(input_path: Path, output_path: Path) -> dict[str, str]:
	"""Strip VSEARCH-style ;size= and ;tax= suffixes from FASTA headers.

	Returns mapping of sanitized_id -> original_header for result annotation.
	"""
	mapping: dict[str, str] = {}
	with input_path.open("r", encoding="utf-8", errors="replace") as src, output_path.open("w", encoding="utf-8") as dst:
		for line in src:
			if line.startswith(">"):
				original = line[1:].strip()
				clean = original.split(";")[0].strip()
				mapping[clean] = original
				dst.write(f">{clean}\n")
			else:
				dst.write(line)
	return mapping


def _remap_bold_output_ids(output_path: Path, header_mapping: dict[str, str]) -> None:
	if not output_path.exists() or not header_mapping:
		return

	def _remap_value(value):
		text = str(value).strip()
		key = text.split(";")[0].strip()
		return header_mapping.get(key, text)

	def _candidate_id_columns(columns: list[str]) -> list[str]:
		preferred = ["ASV_ID", "asv_id", "qseqid", "id", "seq_id", "query"]
		return [column for column in preferred if column in columns]

	# First try structured readers/writers for real parquet/xlsx outputs.
	if pd is not None:
		try:
			if output_path.suffix.lower() == ".parquet":
				df = pd.read_parquet(output_path)
				columns = _candidate_id_columns(list(df.columns))
				if not columns and len(df.columns) > 0:
					columns = [df.columns[0]]
				for column in columns:
					df[column] = df[column].map(_remap_value)
				df.to_parquet(output_path, index=False)
				return
			if output_path.suffix.lower() in {".xlsx", ".xls"}:
				df = pd.read_excel(output_path)
				columns = _candidate_id_columns(list(df.columns))
				if not columns and len(df.columns) > 0:
					columns = [df.columns[0]]
				for column in columns:
					df[column] = df[column].map(_remap_value)
				df.to_excel(output_path, index=False)
				return
		except Exception:
			pass

	# Text fallback for CSV-like test fixtures.
	try:
		text = output_path.read_text(encoding="utf-8")
	except Exception:
		return
	if not text.strip() or "\x00" in text:
		return
	try:
		dialect = csv.Sniffer().sniff(text[:4096]) if text.strip() else csv.excel
	except Exception:
		dialect = csv.excel
	rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
	if not rows:
		return
	fieldnames = list(rows[0].keys())
	columns = _candidate_id_columns(fieldnames)
	if not columns and fieldnames:
		columns = [fieldnames[0]]
	for row in rows:
		for column in columns:
			row[column] = _remap_value(row.get(column, ""))
	with output_path.open("w", newline="", encoding="utf-8") as handle:
		writer = csv.DictWriter(handle, fieldnames=fieldnames)
		writer.writeheader()
		writer.writerows(rows)


def _boldigger_command(fasta_path: Path, db_argument: Path, config: BOLDConfig) -> List[str]:
	if shutil.which("boldigger3"):
		return [
			"boldigger3",
			"identify",
			str(fasta_path),
			str(db_argument),
			"--db",
			str(config.db_source),
			"--mode",
			str(config.mode),
			"--thresholds",
			*[str(value) for value in config.thresholds],
		]
	if BOLDIGGER_PYTHON:
		return [
			BOLDIGGER_PYTHON,
			"-m",
			"boldigger3",
			"identify",
			str(fasta_path),
			str(db_argument),
			"--db",
			str(config.db_source),
			"--mode",
			str(config.mode),
			"--thresholds",
			*[str(value) for value in config.thresholds],
		]
	return [
		"boldigger3",
		"identify",
		str(fasta_path),
		str(db_argument),
		"--db",
		str(config.db_source),
		"--mode",
		str(config.mode),
		"--thresholds",
		*[str(value) for value in config.thresholds],
	]


def run_boldigger(fasta_path: Path, config: BOLDConfig, output_dir: Path) -> StepOutput:
	"""Run BOLDigger3 on a FASTA file and return a StepOutput."""
	output_dir.mkdir(parents=True, exist_ok=True)
	log_path = output_dir / "boldigger.log"
	parquet_path = output_dir / "boldigger_result.parquet"
	xlsx_path = output_dir / "boldigger_result.xlsx"
	sanitized_fasta = output_dir / "sanitized_input.fasta"
	header_mapping = _sanitize_fasta_for_bold(fasta_path, sanitized_fasta)
	log_lines: list[str] = []
	boldigger_available = bool(shutil.which("boldigger3"))
	python_available = bool(BOLDIGGER_PYTHON and (shutil.which(BOLDIGGER_PYTHON) or os.path.isfile(BOLDIGGER_PYTHON)))
	bold_config: BOLDConfig | None = config
	if bold_config is not None:
		ddb_files = list(bold_config.db_path.glob("*.ddb")) if bold_config.db_path.exists() else []
		if not ddb_files and not boldigger_available and not python_available:
			log_lines.append(f"BOLD database empty or missing at {bold_config.db_path} — skipping BOLD.")
			bold_config = None

	if bold_config is None:
		log_path.write_text("\n".join(log_lines), encoding="utf-8")
		return StepOutput(step_name="boldigger", success=False, files={}, log=log_path, error_message=log_lines[-1])

	errors = bold_config.validate()
	if errors:
		log_path.write_text("\n".join(errors), encoding="utf-8")
		return StepOutput(step_name="boldigger", success=False, files={}, log=log_path, error_message="\n".join(errors))

	db_argument = _resolve_db_argument(bold_config)
	command = _boldigger_command(sanitized_fasta, db_argument, bold_config)

	try:
		completed = subprocess.run(command, capture_output=True, text=True, timeout=bold_config.timeout, cwd=str(output_dir), check=True)
		log_path.write_text((completed.stdout or "") + ("\n" if completed.stdout and completed.stderr else "") + (completed.stderr or ""), encoding="utf-8")
	except subprocess.TimeoutExpired:
		log_path.write_text(f"BOLDigger3 timeout after {bold_config.timeout}s", encoding="utf-8")
		return StepOutput(step_name="boldigger", success=False, files={}, log=log_path, error_message=f"BOLDigger3 excedeu timeout de {bold_config.timeout}s.")
	except subprocess.CalledProcessError as exc:
		log_path.write_text((exc.stdout or "") + ("\n" if exc.stdout and exc.stderr else "") + (exc.stderr or ""), encoding="utf-8")
		return StepOutput(step_name="boldigger", success=False, files={}, log=log_path, error_message=f"boldigger3 falhou: {str(exc.stderr or exc.stdout or exc)[:500]}")

	files = {}
	if parquet_path.exists():
		_remap_bold_output_ids(parquet_path, header_mapping)
		files["bold_parquet"] = parquet_path
	if xlsx_path.exists():
		_remap_bold_output_ids(xlsx_path, header_mapping)
		files["bold_xlsx"] = xlsx_path
	if not files:
		return StepOutput(step_name="boldigger", success=False, files={}, log=log_path, error_message="BOLDigger3 executou sem erro mas nao gerou outputs.")
	return StepOutput(step_name="boldigger", success=True, files=files, log=log_path)


def download_bold_db(target_dir: Path) -> bool:
	"""Download the BOLD DuckDB as an explicit setup step."""
	target_dir.mkdir(parents=True, exist_ok=True)
	try:
		subprocess.run(["boldigger3", "download_db", str(target_dir)], check=True, timeout=7200)
		return True
	except Exception:
		return False


def get_bold_reference(marker: str) -> str:
	"""Legacy symbolic helper preserved for compatibility."""
	_ = BOLD_CACHE_ROOT
	return f"bold:{marker.lower()}"