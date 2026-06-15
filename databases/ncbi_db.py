"""NCBI BLAST+ database helpers."""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from core.config import BLASTN_PATH, BLAST_TIMEOUT_SEC, NCBI_BLAST_DB, NCBI_REMOTE_FALLBACK, NCBI_TIMEOUT_SEC
from core.models import StepOutput


@dataclass(frozen=True)
class NCBIConfig:
	local_db_path: Optional[Path] = None
	use_remote: bool = False
	evalue: float = 1e-10
	max_target_seqs: int = 5
	perc_identity: float = 90.0
	qcov_hsp_perc: float = 80.0
	num_threads: int = 8
	word_size: int = 11
	outfmt: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames"
	timeout: int = NCBI_TIMEOUT_SEC

	def validate(self) -> List[str]:
		errors: List[str] = []
		if not shutil.which(BLASTN_PATH) and not os.path.isfile(BLASTN_PATH):
			errors.append("blastn nao encontrado no PATH. Instala BLAST+ ou ajusta BLASTN_PATH.")
		if not self.use_remote and self.local_db_path is None:
			errors.append("NCBIConfig: define local_db_path ou activa use_remote=True.")
		if not self.use_remote and self.local_db_path is not None:
			nal = self.local_db_path.with_suffix(".nal")
			nin = self.local_db_path.with_suffix(".nin")
			if not nal.exists() and not nin.exists():
				errors.append(
					f"Base BLAST local nao encontrada em: {self.local_db_path}\n"
					f"  Esperado: {nal} ou {nin}\n"
					f"  Usa update_blastdb.pl --decompress nt ou define use_remote=True."
				)
		return errors


def _read_fasta_ids(fasta_path: Path) -> set[str]:
	ids: set[str] = set()
	if not fasta_path.exists():
		return ids
	with fasta_path.open("r", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			if line.startswith(">"):
				ids.add(line[1:].split()[0].strip())
	return ids


def _write_no_hit_list(fasta_path: Path, blast_tsv: Path, out_path: Path) -> None:
	fasta_ids = _read_fasta_ids(fasta_path)
	hit_ids: set[str] = set()
	if blast_tsv.exists() and blast_tsv.stat().st_size > 0:
		with blast_tsv.open("r", encoding="utf-8", errors="replace") as handle:
			for line in handle:
				if line.strip() and not line.startswith("#"):
					hit_ids.add(line.split("\t", 1)[0])
	no_hit = sorted(fasta_ids - hit_ids)
	out_path.write_text("\n".join(no_hit), encoding="utf-8")


def run_blast(fasta_path: Path, config: NCBIConfig, output_dir: Path) -> StepOutput:
	"""Run blastn against NCBI nt using local or remote mode."""
	output_dir.mkdir(parents=True, exist_ok=True)
	log_path = output_dir / "blast.log"
	result_path = output_dir / "blast_result.tsv"
	nohit_path = output_dir / "blast_no_hit.txt"

	errors = config.validate()
	# If configured timeout is too low, avoid attempting remote BLAST which
	# would block the pipeline. Require at least 1800s for remote BLAST.
	if config.timeout and config.timeout < 1800:
		return StepOutput(
			step_name="blast",
			success=False,
			files={},
			error_message=(
				"BLAST remoto desactivado: timeout configurado abaixo de 1800s. "
				"Aumenta NCBI_TIMEOUT_SEC em core/config.py ou usa uma base local."
			),
		)
	if errors:
		log_path.write_text("\n".join(errors), encoding="utf-8")
		return StepOutput(step_name="blast", success=False, files={}, log=log_path, error_message="\n".join(errors))

	command = [
		BLASTN_PATH,
		"-query",
		str(fasta_path),
		"-out",
		str(result_path),
		"-outfmt",
		config.outfmt,
		"-evalue",
		str(config.evalue),
		"-max_target_seqs",
		str(config.max_target_seqs),
		"-perc_identity",
		str(config.perc_identity),
		"-qcov_hsp_perc",
		str(config.qcov_hsp_perc),
		"-word_size",
		str(config.word_size),
	]

	if config.use_remote:
		command.extend(["-db", "nt", "-remote"])
	else:
		command.extend(["-db", str(config.local_db_path), "-num_threads", str(config.num_threads)])

	try:
		completed = subprocess.run(command, capture_output=True, text=True, check=True, timeout=config.timeout)
		log_path.write_text((completed.stdout or "") + ("\n" if completed.stdout and completed.stderr else "") + (completed.stderr or ""), encoding="utf-8")
	except subprocess.TimeoutExpired:
		log_path.write_text(f"BLAST timeout after {config.timeout}s", encoding="utf-8")
		return StepOutput(step_name="blast", success=False, files={}, log=log_path, error_message=f"BLAST excedeu o timeout de {config.timeout}s.")
	except subprocess.CalledProcessError as exc:
		log_path.write_text((exc.stdout or "") + ("\n" if exc.stdout and exc.stderr else "") + (exc.stderr or ""), encoding="utf-8")
		return StepOutput(step_name="blast", success=False, files={}, log=log_path, error_message=f"blastn falhou: {str(exc.stderr or exc.stdout or exc)[:500]}")

	if not result_path.exists():
		result_path.write_text("", encoding="utf-8")
	_write_no_hit_list(fasta_path, result_path, nohit_path)
	hit_count = 0.0
	if result_path.exists() and result_path.stat().st_size > 0:
		with result_path.open("r", encoding="utf-8", errors="replace") as handle:
			for line in handle:
				if line.strip() and not line.startswith("#"):
					hit_count += 1.0
	return StepOutput(
		step_name="blast",
		success=True,
		files={"blast_tsv": result_path, "blast_no_hit": nohit_path, "log": log_path},
		metrics={"hit_count": hit_count},
	)


def download_blast_db(target_dir: Path, db_name: str = "nt") -> bool:
	"""Download the NCBI BLAST DB with update_blastdb.pl as an explicit setup step."""
	if shutil.which("update_blastdb.pl") is None:
		return False
	target_dir.mkdir(parents=True, exist_ok=True)
	try:
		subprocess.run(["update_blastdb.pl", "--decompress", db_name], cwd=str(target_dir), check=True, timeout=86400)
		return True
	except Exception:
		return False


def get_ncbi_reference(marker: str) -> str:
	"""Legacy symbolic helper preserved for compatibility."""
	_ = (NCBI_BLAST_DB, NCBI_REMOTE_FALLBACK)
	return f"ncbi:{marker.lower()}"