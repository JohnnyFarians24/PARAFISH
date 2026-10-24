"""NCBI BLAST+ database helpers."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import threading
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


# How a remote search reaches NCBI. "urlapi" sends the sequences through the BLAST URL API
# (Biopython's qblast), the service behind the BLAST web page; "blastn" runs
# ``blastn -remote``. Both are documented NCBI routes, but on the evaluation laptop a
# single 388 bp OTU sent with ``blastn -remote`` had no answer after 15 min against nt and
# after 5 min against core_nt, while the same OTU came back in about 1 min through the URL
# API, and in seconds on the web page.
NCBI_REMOTE_METHOD = (os.getenv("NCBI_REMOTE_METHOD", "urlapi").strip().lower() or "urlapi")
# core_nt is the default nucleotide database of web BLAST: GenBank without most of the
# whole-genome and unplaced sequence that makes nt more than twice its size.
NCBI_REMOTE_DB = (os.getenv("NCBI_REMOTE_DB", "core_nt").strip() or "core_nt")
# Sequences per URL API request. One request per run is what NCBI asks for, but a request
# with thousands of OTUs is slow to queue and loses everything if it fails.
try:
	NCBI_REMOTE_BATCH = max(1, int(os.getenv("NCBI_REMOTE_BATCH", "100") or 100))
except ValueError:
	NCBI_REMOTE_BATCH = 100
# NCBI: "Do not contact the server more often than once every 10 seconds."
REMOTE_REQUEST_GAP_SEC = 10


def blast_db_exists(path) -> bool:
	"""True when ``path`` names a BLAST database: the prefix itself, or its .nal/.nin.

	A BLAST database is a prefix (``SSU_eukaryote_rRNA``) followed by several files
	(``.nin``, ``.nsq``, ...); there is no file with the bare name, so ``Path.exists()``
	alone reported every real database as missing and the run fell back to remote.
	"""
	if not path:
		return False
	p = Path(path)
	return p.exists() or p.with_suffix(".nal").exists() or p.with_suffix(".nin").exists() or Path(str(p) + ".nal").exists() or Path(str(p) + ".nin").exists()


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
	# Remote route and database (see NCBI_REMOTE_METHOD and NCBI_REMOTE_DB above).
	remote_method: str = NCBI_REMOTE_METHOD
	remote_db: str = NCBI_REMOTE_DB

	def effective_timeout(self) -> int:
		"""Seconds to allow this BLAST. Local searches never inherit the remote gate."""
		if self.timeout and self.timeout > 0:
			return self.timeout
		return BLAST_TIMEOUT_SEC

	def validate(self) -> List[str]:
		errors: List[str] = []
		needs_blastn = not (self.use_remote and self.remote_method == "urlapi")
		if needs_blastn and not shutil.which(BLASTN_PATH) and not os.path.isfile(BLASTN_PATH):
			errors.append("blastn not found on PATH. Install BLAST+ or adjust BLASTN_PATH.")
		if not self.use_remote and self.local_db_path is None:
			errors.append("NCBIConfig: set local_db_path or enable use_remote=True.")
		if not self.use_remote and self.local_db_path is not None:
			nal = self.local_db_path.with_suffix(".nal")
			nin = self.local_db_path.with_suffix(".nin")
			if not nal.exists() and not nin.exists():
				errors.append(
					f"Local BLAST database not found at: {self.local_db_path}\n"
					f"  Expected: {nal} or {nin}\n"
					f"  Use update_blastdb.pl --decompress nt, or set use_remote=True."
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


def _read_fasta(fasta_path: Path) -> list[tuple[str, str]]:
	"""``[(header, sequence), ...]`` in file order."""
	records: list[tuple[str, str]] = []
	header, chunks = None, []
	with fasta_path.open("r", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			line = line.rstrip("\n")
			if line.startswith(">"):
				if header is not None:
					records.append((header, "".join(chunks)))
				header, chunks = line[1:].strip(), []
			elif line.strip():
				chunks.append(line.strip())
	if header is not None:
		records.append((header, "".join(chunks)))
	return records


def _call_with_timeout(func, timeout: float):
	"""``(result, exception, timed_out)`` of ``func()`` run in a worker thread.

	qblast polls NCBI until the search is ready and has no time limit of its own; the
	worker is a daemon, so a search that never returns cannot keep the app alive.
	"""
	box: dict = {}

	def target():
		try:
			box["result"] = func()
		except Exception as exc:  # reported to the caller, never swallowed
			box["error"] = exc

	worker = threading.Thread(target=target, daemon=True)
	worker.start()
	worker.join(timeout)
	if worker.is_alive():
		return None, None, True
	return box.get("result"), box.get("error"), False


def _qblast_xml(sequences: str, config: NCBIConfig) -> str:
	"""One URL API search; returns the BLAST XML text. Isolated so tests can replace it."""
	from Bio.Blast import NCBIWWW
	from core.config import get_entrez_config

	email = get_entrez_config().get("email") or ""
	if email:
		NCBIWWW.email = email
	kwargs = dict(
		megablast=(str(config.task) == "megablast"),
		hitlist_size=int(config.max_target_seqs),
		expect=float(config.evalue),
		format_type="XML",
	)
	if config.word_size:
		kwargs["word_size"] = int(config.word_size)
	handle = NCBIWWW.qblast("blastn", str(config.remote_db), sequences, **kwargs)
	try:
		return handle.read()
	finally:
		handle.close()


def _accession_taxids(accessions: list[str]) -> dict[str, str]:
	"""``{accession.version: taxid}`` from Entrez esummary. Best-effort.

	The XML of the URL API names the hit but carries no taxid, and the taxid is what the
	lineage lookup (ncbi_lineage) and the consensus need.
	"""
	found: dict[str, str] = {}
	if not accessions:
		return found
	try:
		from Bio import Entrez
		from core.config import get_entrez_config
	except ImportError:
		return found
	cfg = get_entrez_config()
	if not cfg.get("email"):
		return found
	Entrez.email = cfg["email"]
	Entrez.tool = "PARAFISH"
	if cfg.get("api_key"):
		Entrez.api_key = cfg["api_key"]
	unique = list(dict.fromkeys(accessions))
	for start in range(0, len(unique), 200):
		chunk = unique[start : start + 200]
		try:
			handle = Entrez.esummary(db="nuccore", id=",".join(chunk))
			records = Entrez.read(handle)
			handle.close()
		except Exception:
			continue
		for rec in records:
			acc = str(rec.get("AccessionVersion") or rec.get("Caption") or "").strip()
			# Entrez.read gives IntegerElement objects, whose str() is their repr
			# ("IntegerElement(8063, attributes={})"), not the number.
			try:
				taxid = str(int(rec.get("TaxId") or 0))
			except (TypeError, ValueError):
				taxid = ""
			if acc and taxid and taxid != "0":
				found[acc] = taxid
				found.setdefault(acc.split(".")[0], taxid)
	return found


def _taxid_names(taxids: list[str]) -> dict[str, str]:
	"""``{taxid: species name}`` through the cached lineage lookup. Best-effort."""
	if not taxids:
		return {}
	try:
		from databases.ncbi_lineage import resolve_lineages
	except ImportError:
		return {}
	lineages, _warnings = resolve_lineages(taxids)
	return {tid: (ranks.get("species") or ranks.get("genus") or "") for tid, ranks in lineages.items()}


def _gap_opens(aligned: str) -> int:
	return len(re.findall(r"-+", aligned or ""))


def _xml_to_rows(xml_text: str, config: NCBIConfig) -> list[dict]:
	"""Hit rows from BLAST XML, one per HSP kept, with blastn's tabular field names.

	blastn applies -perc_identity and -qcov_hsp_perc itself; the URL API does not, so the
	same thresholds are applied here and the two routes keep the same evidence standard.
	"""
	from io import StringIO

	from Bio.Blast import NCBIXML

	rows: list[dict] = []
	for record in NCBIXML.parse(StringIO(xml_text)):
		qseqid = (record.query or record.query_id or "").split()[0]
		qlen = int(record.query_length or record.query_letters or 0)
		for alignment in record.alignments[: int(config.max_target_seqs)]:
			hit_id = re.sub(r"^gi\|\d+\|", "", alignment.hit_id or "")
			for hsp in alignment.hsps:
				length = int(hsp.align_length or 0)
				if not length:
					continue
				ident = int(hsp.identities or 0)
				gaps = int(hsp.gaps) if hsp.gaps is not None else (hsp.query or "").count("-") + (hsp.sbjct or "").count("-")
				pident = 100.0 * ident / length
				qcov = 100.0 * (abs(int(hsp.query_end) - int(hsp.query_start)) + 1) / qlen if qlen else 0.0
				if pident < float(config.perc_identity) or qcov < float(config.qcov_hsp_perc):
					continue
				rows.append(
					{
						"qseqid": qseqid,
						"sseqid": hit_id,
						"accession": str(alignment.accession or "").strip(),
						"pident": f"{pident:.3f}",
						"qcovhsp": str(int(round(qcov))),
						"length": str(length),
						"mismatch": str(max(0, length - ident - gaps)),
						"gapopen": str(_gap_opens(hsp.query) + _gap_opens(hsp.sbjct)),
						"qstart": str(hsp.query_start),
						"qend": str(hsp.query_end),
						"sstart": str(hsp.sbjct_start),
						"send": str(hsp.sbjct_end),
						"evalue": f"{float(hsp.expect):.3g}",
						"bitscore": f"{float(hsp.bits):.1f}",
						"staxids": "",
						# Left blank, not "N/A": a blank species is filled from the taxid lineage,
						# while "N/A" would be taken as a species name.
						"sscinames": "",
						"scomnames": "",
					}
				)
	return rows


def _run_remote_urlapi(fasta_path: Path, config: NCBIConfig, result_path: Path, nohit_path: Path, log_path: Path) -> StepOutput:
	"""Remote search through the BLAST URL API, written in blastn's tabular format.

	The output file has exactly the columns of ``config.outfmt``, so every reader of
	``blast_result.tsv`` downstream works unchanged whichever route produced it.
	"""
	log_lines = [f"NCBI URL API (qblast): db={config.remote_db}, task={config.task}, batch={NCBI_REMOTE_BATCH}"]
	try:
		import Bio.Blast.NCBIWWW  # noqa: F401
	except ImportError:
		log_path.write_text("Biopython is not installed.", encoding="utf-8")
		return StepOutput(step_name="blast", success=False, files={}, log=log_path, error_message="Biopython is required for the NCBI URL API route (NCBI_REMOTE_METHOD=urlapi).")

	records = _read_fasta(fasta_path)
	fields = config.outfmt.split()[1:]
	deadline = time.monotonic() + config.effective_timeout()
	rows: list[dict] = []
	for start in range(0, len(records), NCBI_REMOTE_BATCH):
		batch = records[start : start + NCBI_REMOTE_BATCH]
		if start:
			time.sleep(REMOTE_REQUEST_GAP_SEC)
		remaining = deadline - time.monotonic()
		if remaining <= 0:
			remaining = 0
		fasta_text = "".join(f">{header}\n{seq}\n" for header, seq in batch)
		xml_text, error, timed_out = _call_with_timeout(lambda: _qblast_xml(fasta_text, config), remaining)
		label = f"batch {start // NCBI_REMOTE_BATCH + 1} ({len(batch)} sequences)"
		if timed_out:
			log_lines.append(f"{label}: no answer before the time limit")
			log_path.write_text("\n".join(log_lines), encoding="utf-8")
			return StepOutput(
				step_name="blast",
				success=False,
				files={},
				log=log_path,
				error_message=(
					f"NCBI did not answer within {config.effective_timeout()}s (URL API, db={config.remote_db}). "
					f"Raise NCBI_TIMEOUT_SEC in .env, try again later, or use a local BLAST database."
				),
			)
		if error is not None:
			text = str(error)
			log_lines.append(f"{label}: {text[:400]}")
			log_path.write_text("\n".join(log_lines), encoding="utf-8")
			refusal = _server_refusal(text)
			detail = (
				f"NCBI refused the search (URL API): {refusal}" if refusal else f"NCBI URL API search failed: {text[:400]}"
			)
			return StepOutput(step_name="blast", success=False, files={}, log=log_path, error_message=detail)
		batch_rows = _xml_to_rows(xml_text or "", config)
		log_lines.append(f"{label}: {len(batch_rows)} hits")
		rows.extend(batch_rows)

	taxids = _accession_taxids([row["accession"] for row in rows if row["accession"]])
	names = _taxid_names(sorted(set(taxids.values())))
	for row in rows:
		tid = taxids.get(row["accession"]) or taxids.get(row["accession"].split(".")[0], "")
		row["staxids"] = tid
		if names.get(tid):
			row["sscinames"] = names[tid]
	missing = sum(1 for row in rows if not row["staxids"])
	if missing:
		log_lines.append(f"{missing} hit(s) without a taxid (Entrez unavailable or ENTREZ_EMAIL not set)")

	with result_path.open("w", encoding="utf-8") as handle:
		for row in rows:
			handle.write("\t".join(str(row.get(field, "")) for field in fields) + "\n")
	log_path.write_text("\n".join(log_lines), encoding="utf-8")
	_write_no_hit_list(fasta_path, result_path, nohit_path)
	return StepOutput(
		step_name="blast",
		success=True,
		files={"blast_tsv": result_path, "blast_no_hit": nohit_path, "log": log_path},
		metrics={"hit_count": float(len(rows))},
	)


def run_blast(fasta_path: Path, config: NCBIConfig, output_dir: Path) -> StepOutput:
	"""Run a BLAST search against NCBI: local database, URL API (default remote) or blastn -remote."""
	output_dir.mkdir(parents=True, exist_ok=True)
	log_path = output_dir / "blast.log"
	result_path = output_dir / "blast_result.tsv"
	nohit_path = output_dir / "blast_no_hit.txt"

	errors = config.validate()
	if errors:
		log_path.write_text("\n".join(errors), encoding="utf-8")
		return StepOutput(step_name="blast", success=False, files={}, log=log_path, error_message="\n".join(errors))

	if config.use_remote and config.remote_method == "urlapi":
		return _run_remote_urlapi(fasta_path, config, result_path, nohit_path, log_path)

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
		command.extend(["-db", str(config.remote_db), "-remote"])
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
					f"BLAST exceeded the {timeout}s timeout. Remote BLAST runs on NCBI's servers "
					f"and takes about 1 min per sequence, plus queueing. Raise NCBI_TIMEOUT_SEC "
					f"in .env, or reduce the number of sequences."
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
				"the installed blastn cannot talk to NCBI (TLS/MBEDTLS error, typical of the "
				"Ubuntu blast 2.9.0 package). Install the official BLAST+ and point BLASTN_PATH at it: "
				"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
			)
		elif "cpu usage limit" in last_detail.lower() or "sigxcpu" in last_detail.lower():
			detail = (
				f"NCBI refused the search for exceeding the public CPU quota (SIGXCPU) across "
				f"{len(REMOTE_BACKOFF_SEC)} spaced-out attempts. This is a server REFUSAL, not "
				f"an absence of hits. task='{config.task}'"
				+ (f", word_size={config.word_size}" if config.word_size else "")
				+ ". More sensitive searches (task=blastn/dc-megablast, low word_size) against the whole "
				"nt database are refused by the public quota; use task=megablast, or a local BLAST database."
			)
		else:
			detail = f"blastn failed: {last_detail[:400]}"
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