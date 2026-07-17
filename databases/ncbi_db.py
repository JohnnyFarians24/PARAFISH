"""NCBI BLAST+ database helpers."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from core.process import find_failure_signature
from core.config import BLASTDB_PATH, BLASTN_PATH, BLAST_TIMEOUT_SEC, NCBI_BLAST_DB, NCBI_REMOTE_FALLBACK, NCBI_TIMEOUT_SEC
from core.models import StepOutput


# NCBI kills a remote search that costs too much CPU on its shared public quota, and it
# does so per query, with this in stderr. Measured on the reference run's 9 OTUs against
# nt (all three tried, same query, same day):
#
#   -task blastn -word_size 11   -> REFUSED, SIGXCPU on all 9        (what we used to send)
#   -task dc-megablast           -> REFUSED, SIGXCPU on all 9
#   -task megablast (default)    -> 123 hits, 9/9 OTUs, ~3 min       <- the only one that runs
#
# So the throttling was not arbitrary: word_size 11 against all of nt is a full blastn
# search, the most expensive thing you can ask for, and the public quota exists to stop
# exactly that. This is OUR request to fix, not NCBI's mood to wait out.
BLAST_FAILURE_SIGNATURES = (
	"cpu usage limit",
	"sigxcpu",
	"search failed",
	"failed to receive reply",
	"mbedtls",
	"cannot connect",
	"service unavailable",
)


# Waits between retries of a REFUSED remote search, in seconds. Spaced and growing: the
# public quota is shared, so load varies and a later attempt can succeed — but the first
# wait already clears NCBI's "no more than one request per 10 seconds", and after the last
# one we stop and say so. Never a silent zero.
REMOTE_BACKOFF_SEC = (30, 120, 0)


def _server_refusal(text: str) -> Optional[str]:
	"""The first line of ``text`` that shows NCBI refused the search, or ``None``.

	blastn exits 0 having printed these and written an empty result file, so nothing
	downstream can tell a refusal from "this sample contains nothing in nt" unless it
	reads the log. That is how a killed search reached source_status as ``no_hits``.
	"""
	return find_failure_signature(text, BLAST_FAILURE_SIGNATURES)


@dataclass(frozen=True)
class NCBIConfig:
	local_db_path: Optional[Path] = None
	use_remote: bool = False
	evalue: float = 1e-10
	max_target_seqs: int = 5
	perc_identity: float = 90.0
	qcov_hsp_perc: float = 80.0
	num_threads: int = 8
	# BLAST task. Remote runs on NCBI's shared quota and only megablast completes there
	# (see above). A LOCAL database is our own CPU, so it can afford the sensitive search
	# that finds a parasite at ~90% identity, which megablast's 28 bp exact seed can miss.
	task: str = "megablast"
	# None = let the task pick its own seed length (28 for megablast, 11 for blastn).
	# Pinning 11 while leaving the task at megablast is what produced a search with
	# blastn's cost and no way to complete.
	word_size: Optional[int] = None
	# ``qcovhsp`` is what lets the UI report query coverage next to identity: 99%
	# identity over 60 bp of a 400 bp amplicon is not the evidence that 99% over the
	# whole query is, and identity alone cannot tell the two apart.
	outfmt: str = "6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames"
	# 0 means "not configured", not "give up immediately". NCBI_TIMEOUT_SEC defaults
	# to 0 to keep *remote* BLAST opt-in, and that zero used to be handed straight to
	# subprocess(timeout=0), which killed every local BLAST on the spot.
	timeout: int = NCBI_TIMEOUT_SEC

	def effective_timeout(self) -> int:
		"""Seconds to allow this BLAST. Local searches never inherit the remote gate."""
		if self.timeout and self.timeout > 0:
			return self.timeout
		return BLAST_TIMEOUT_SEC

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
		"-task",
		str(config.task),
		"-evalue",
		str(config.evalue),
		"-max_target_seqs",
		str(config.max_target_seqs),
		"-perc_identity",
		str(config.perc_identity),
		"-qcov_hsp_perc",
		str(config.qcov_hsp_perc),
	]
	if config.word_size:
		command.extend(["-word_size", str(config.word_size)])

	if config.use_remote:
		# -num_threads is meaningless (and warned about) for a remote search: the
		# search runs on NCBI's hardware.
		#
		# Every OTU goes in ONE request: -query is the whole multi-FASTA and blastn
		# submits it as a single job. Per-sequence requests would be both slower and the
		# kind of traffic NCBI's limits exist to stop.
		command.extend(["-db", "nt", "-remote"])
	else:
		command.extend(["-db", str(config.local_db_path), "-num_threads", str(config.num_threads)])

	# blastn resolves taxids to species names through the taxdb files, which it finds
	# via $BLASTDB. Without it every `sscinames` is "N/A" — the alignment is fine but
	# carries no taxonomy, which is worse than useless for a classification step.
	env = os.environ.copy()
	if BLASTDB_PATH:
		env["BLASTDB"] = BLASTDB_PATH

	# NCBI asks callers not to contact the server more than once every 10 s. One run makes
	# one request, so the only way we approach that is by retrying — hence spaced, growing
	# waits rather than an immediate retry. A refused search is retried because the public
	# quota is shared and load varies; it is NOT retried forever, and it never degrades
	# into a silent zero.
	attempts_log: list[str] = []
	last_detail = ""
	for attempt, backoff in enumerate(REMOTE_BACKOFF_SEC if config.use_remote else (0,), start=1):
		try:
			timeout = config.effective_timeout()
			completed = subprocess.run(command, capture_output=True, text=True, check=True, timeout=timeout, env=env)
			output = (completed.stdout or "") + ("\n" if completed.stdout and completed.stderr else "") + (completed.stderr or "")
		except subprocess.TimeoutExpired:
			timeout = config.effective_timeout()
			log_path.write_text("\n".join([*attempts_log, f"BLAST timeout after {timeout}s"]), encoding="utf-8")
			return StepOutput(
				step_name="blast",
				success=False,
				files={},
				log=log_path,
				error_message=(
					f"BLAST excedeu o timeout de {timeout}s. O BLAST remoto corre nos servidores do NCBI "
					f"e demora cerca de 1 min por sequência, mais fila de espera. Aumenta NCBI_TIMEOUT_SEC "
					f"no .env, ou reduz o número de sequências."
				),
			)
		except subprocess.CalledProcessError as exc:
			output = (exc.stdout or "") + ("\n" if exc.stdout and exc.stderr else "") + (exc.stderr or "")
			refusal = _server_refusal(output) or str(exc)
			last_detail = refusal
			attempts_log.append(f"--- attempt {attempt} (exit {exc.returncode}) ---\n{output}")
			break

		# blastn EXITS 0 having printed "Severe Error ... SIGXCPU" for every query and
		# written an empty result. Without this the step returns success with hit_count=0
		# and the refusal reaches source_status wearing the word "no_hits".
		refusal = _server_refusal(output)
		attempts_log.append(f"--- attempt {attempt} ---\n{output}")
		if refusal is None:
			log_path.write_text("\n".join(attempts_log), encoding="utf-8")
			break
		last_detail = refusal
		if backoff and attempt < len(REMOTE_BACKOFF_SEC):
			attempts_log.append(f"NCBI refused the search; waiting {backoff}s before retry {attempt + 1}.")
			time.sleep(backoff)
	else:
		last_detail = last_detail or "NCBI refused the search"

	if last_detail:
		log_path.write_text("\n".join(attempts_log), encoding="utf-8")
		# The Ubuntu blastn 2.9.0 package ships a broken TLS stack and can never reach
		# NCBI; no timeout or retry helps, so say what actually needs to change.
		if "MBEDTLS" in last_detail or "Failed to receive reply" in last_detail:
			detail = (
				"o blastn instalado nao consegue falar com o NCBI (erro de TLS/MBEDTLS, tipico do "
				"pacote blast 2.9.0 do Ubuntu). Instala o BLAST+ oficial e aponta BLASTN_PATH para ele: "
				"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
			)
		elif "cpu usage limit" in last_detail.lower() or "sigxcpu" in last_detail.lower():
			detail = (
				f"o NCBI recusou a pesquisa por exceder a quota de CPU pública (SIGXCPU) em "
				f"{len(REMOTE_BACKOFF_SEC)} tentativas espaçadas. Isto é uma RECUSA do servidor, não "
				f"ausência de hits. task='{config.task}'"
				+ (f", word_size={config.word_size}" if config.word_size else "")
				+ ". Buscas mais sensíveis (task=blastn/dc-megablast, word_size baixo) contra o nt "
				"inteiro são recusadas pela quota pública; usa task=megablast, ou uma base BLAST local."
			)
		else:
			detail = f"blastn falhou: {last_detail[:400]}"
		return StepOutput(step_name="blast", success=False, files={}, log=log_path, error_message=detail)

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