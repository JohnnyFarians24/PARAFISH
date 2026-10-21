"""BOLD reference database helpers."""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import List, Optional, Tuple

try:
	import pandas as pd
except Exception:  # pragma: no cover
	pd = None

from core.config import BOLD_CACHE_ROOT, BOLDIGGER_PYTHON, BOLDIGGER_TIMEOUT_SEC, get_bold_credentials
from core.models import StepOutput


# boldigger3 < 2.0 crashes with ``TypeError: 'int' object is not callable`` in
# id_engine.main() (it calls ``pbar.total(...)`` instead of assigning to it)
# whenever BOLD requeues an incomplete download — which happens on essentially
# every real run. The bug is present in every 1.x release up to 1.6.3 and was
# fixed upstream in the 2.x/3.x line, so 2.0 is the hard minimum.
MIN_BOLDIGGER_VERSION: Tuple[int, ...] = (2, 0)


def _boldigger_base_command() -> List[str]:
	# An explicitly configured interpreter wins over PATH discovery: BOLDIGGER_PYTHON
	# exists precisely to point at a separate env, which it could never do while a
	# stale ``boldigger3`` on PATH took priority.
	if BOLDIGGER_PYTHON:
		return [BOLDIGGER_PYTHON, "-m", "boldigger3"]
	if shutil.which("boldigger3"):
		return ["boldigger3"]
	return ["boldigger3"]


def boldigger_available() -> bool:
	return bool(shutil.which("boldigger3")) or bool(
		BOLDIGGER_PYTHON and (shutil.which(BOLDIGGER_PYTHON) or os.path.isfile(BOLDIGGER_PYTHON))
	)


@lru_cache(maxsize=1)
def boldigger_version() -> Optional[Tuple[int, ...]]:
	"""Return the installed boldigger3 version, or ``None`` if undetectable."""
	if not boldigger_available():
		return None
	try:
		completed = subprocess.run(
			[*_boldigger_base_command(), "--version"],
			capture_output=True,
			text=True,
			timeout=120,
		)
	except Exception:
		# An undetectable version must never block a run; callers treat None as
		# "assume the installed CLI is fine".
		return None
	streams = "\n".join(part for part in (completed.stdout, completed.stderr) if isinstance(part, str))
	# boldigger3 prints the bare version on a line of its own, but surrounds it with
	# noise (dependency-install chatter, "version is outdated", tracebacks). Only a
	# whole-line match is safe: a loose search happily picks the "3.12" out of a
	# ``.../lib/python3.12/site-packages/...`` path in a traceback.
	match = None
	for line in streams.splitlines():
		found = re.fullmatch(r"\s*(?:boldigger3\s+)?(\d+)\.(\d+)(?:\.(\d+))?\s*", line)
		if found:
			match = found
	if match is None:
		return None
	return tuple(int(part) for part in match.groups() if part is not None)


def _format_version(version: Optional[Tuple[int, ...]]) -> str:
	return ".".join(str(part) for part in version) if version else "desconhecida"


def find_local_bold_db(db_path: Path) -> Optional[Path]:
	"""Return the local BOLD DuckDB (``*.ddb``) file, if one has been downloaded."""
	try:
		if db_path.is_file() and db_path.suffix == ".ddb":
			return db_path
		if db_path.is_dir():
			return next(iter(sorted(db_path.glob("*.ddb"))), None)
	except OSError:
		return None
	return None


# Exactly the columns boldigger3's ``download_db`` projects out of the BOLD parquet.
# add_metadata LEFT JOINs against them, so the table must exist with these names even
# when we have no rows to put in it.
_BOLD_METADATA_COLUMNS = (
	"processid VARCHAR",
	"sex VARCHAR",
	"life_stage VARCHAR",
	"inst VARCHAR",
	'"country/ocean" VARCHAR',
	"identified_by VARCHAR",
	"identification_method VARCHAR",
	"coord VARCHAR",
	"nuc VARCHAR",
	"marker_code VARCHAR",
)


def bold_metadata_db_status(ddb_path: Path) -> tuple[bool, int]:
	"""``(has_bold_public_table, row_count)`` for a local BOLD DuckDB."""
	try:
		import duckdb
	except ImportError:
		return False, 0
	try:
		with duckdb.connect(str(ddb_path), read_only=True) as con:
			tables = [name for (name,) in con.execute("SHOW TABLES").fetchall()]
			if "bold_public" not in tables:
				return False, 0
			return True, int(con.execute("SELECT count(*) FROM bold_public").fetchone()[0])
	except Exception:
		return False, 0


def ensure_bold_metadata_db(db_dir: Path) -> tuple[Optional[Path], Optional[str]]:
	"""Guarantee a ``*.ddb`` that boldigger3 v3 can actually open. Returns (path, note).

	boldigger3 v3 splits the work: ``id_engine`` fetches the hits (taxonomy and
	% identity) from BOLD's remote API, and only ``add_metadata`` touches this local
	database — purely to LEFT JOIN cosmetic specimen fields (sex, institution,
	country) onto them.

	Building the real database means materialising a 2.3 GB parquet into DuckDB,
	which the OOM killer terminates on a small machine — leaving behind the empty,
	table-less ``.ddb`` we found here. ``add_metadata`` then aborts on
	"table bold_public does not exist", destroying an identification that had already
	succeeded.

	An empty-but-well-formed ``bold_public`` costs nothing and makes that LEFT JOIN
	resolve to NULLs, so the run completes with full taxonomy and identity. The only
	thing lost is metadata the pipeline never consumed.
	"""
	try:
		import duckdb
	except ImportError:
		return find_local_bold_db(db_dir), None

	db_dir = Path(db_dir)
	existing = find_local_bold_db(db_dir)
	if existing is not None:
		has_table, rows = bold_metadata_db_status(existing)
		if has_table and rows > 0:
			return existing, None
		if has_table:
			return existing, "Local BOLD database has no records: specimen metadata will come back empty."
		target = existing
		note = (
			f"The local BOLD database ({existing.name}) is empty — the BOLD download finished the "
			f"parquet but ran out of memory building the DuckDB. An empty metadata table was "
			f"created so identification can still run: the taxonomy and the % identity "
			f"come from the remote BOLD API and are unaffected; only the specimen metadata "
			f"(sex, institution, country) is left blank."
		)
	else:
		db_dir.mkdir(parents=True, exist_ok=True)
		target = db_dir / "BOLD_metadata_empty.ddb"
		note = (
			"No local BOLD database. An empty metadata table was created: identification "
			"(taxonomy and % identity) comes from the remote BOLD API and works regardless."
		)

	try:
		with duckdb.connect(str(target)) as con:
			con.execute(f"CREATE TABLE IF NOT EXISTS bold_public ({', '.join(_BOLD_METADATA_COLUMNS)})")
	except Exception as exc:
		return existing, f"Could not prepare the local BOLD database: {exc}"
	return target, note


@dataclass(frozen=True)
class BOLDConfig:
	db_path: Path
	db_source: int = 1
	mode: int = 2
	thresholds: List[int] = field(default_factory=lambda: [97, 95, 90, 85])
	timeout: int = BOLDIGGER_TIMEOUT_SEC

	def validate(self) -> List[str]:
		errors: List[str] = []
		if not boldigger_available():
			errors.append("boldigger3 not found on PATH and BOLDIGGER_PYTHON is not set. Install with: pip install 'boldigger3>=3'")
		else:
			# When the version cannot be determined, assume the caller knows what it
			# is doing rather than blocking the run on a guess.
			version = boldigger_version()
			if version is not None and version < MIN_BOLDIGGER_VERSION:
				errors.append(
					f"boldigger3 {_format_version(version)} has a known bug (id_engine calls 'pbar.total(...)') "
					f"that crashes identification as soon as BOLD re-queues a download. "
					f"Update with: pip install -U 'boldigger3>=3'"
				)
		# NOTE: deliberately no local-*.ddb check. ``identify`` queries the remote
		# BOLD v5 engine and never opens a local database, so requiring one here only
		# blocked runs that would have worked.
		if self.db_source not in range(1, 9):
			errors.append(f"db_source deve ser 1-8, recebido: {self.db_source}")
		if self.mode not in (1, 2, 3):
			errors.append(f"mode deve ser 1, 2 ou 3, recebido: {self.mode}")
		if len(self.thresholds) != 4:
			errors.append(f"thresholds deve ter 4 valores, recebido: {self.thresholds}")
		return errors


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

	name = output_path.name.lower()
	is_parquet = name.endswith(".parquet") or name.endswith(".parquet.snappy")

	# First try structured readers/writers for real parquet/xlsx outputs.
	if pd is not None:
		try:
			if is_parquet:
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


def _boldigger_command(fasta_path: Path, config: BOLDConfig) -> List[str]:
	"""Build the ``boldigger3 identify`` command line for the installed version.

	The two CLI generations are not compatible, and which one you get depends on
	which interpreter ``BOLDIGGER_PYTHON`` points at:

	* 3.x — ``identify <fasta> <db_path>``. Hits come from the remote BOLD v5 engine,
	  but specimen metadata is joined from a locally downloaded DuckDB, so the
	  ``*.ddb`` positional is REQUIRED and argparse fails without it.
	* 1.x — ``identify <fasta>``. Fully remote; an extra positional is rejected with
	  "unrecognized arguments".
	"""
	positionals = [str(fasta_path)]
	version = boldigger_version()
	if version is not None and version >= MIN_BOLDIGGER_VERSION:
		local_db = find_local_bold_db(config.db_path)
		if local_db is not None:
			positionals.append(str(local_db))
	return [
		*_boldigger_base_command(),
		"identify",
		*positionals,
		"--db",
		str(config.db_source),
		"--mode",
		str(config.mode),
		"--thresholds",
		*[str(value) for value in config.thresholds],
	]


def _find_boldigger_outputs(output_dir: Path, fasta_stem: str) -> dict[str, Path]:
	"""Locate the files boldigger3 actually wrote.

	boldigger3 3.x writes into a ``boldigger3_data`` subdirectory next to the
	input FASTA, while 1.x wrote straight beside it. Both layouts (and the old
	``boldigger_result.*`` names) are probed so the runner survives an upgrade.
	"""
	data_dir = output_dir / "boldigger3_data"
	found: dict[str, Path] = {}
	candidates = {
		"bold_parquet": [
			data_dir / f"{fasta_stem}_identification_result.parquet.snappy",
			output_dir / f"{fasta_stem}_identification_result.parquet.snappy",
			output_dir / f"{fasta_stem}_identification_result.parquet",
			output_dir / "boldigger_result.parquet",
		],
		"bold_xlsx": [
			data_dir / f"{fasta_stem}_identification_result.xlsx",
			output_dir / f"{fasta_stem}_identification_result.xlsx",
			output_dir / "boldigger_result.xlsx",
		],
	}
	for key, paths in candidates.items():
		for path in paths:
			if path.exists() and path.stat().st_size > 0:
				found[key] = path
				break
	if not found:
		# Last resort: any identification result the tool may have named differently.
		for path in sorted(output_dir.rglob("*identification_result*")):
			if path.suffix.lower() in {".xlsx", ".xls"}:
				found.setdefault("bold_xlsx", path)
			elif "parquet" in path.name.lower():
				found.setdefault("bold_parquet", path)
	return found


# BOLD's public portal answers 503 whenever it is busy, and boldigger3's
# additional-data phase maps only {200, 404, 500} — so a 503 raises KeyError and
# takes the whole run down. The identification itself has already succeeded and been
# persisted by then, so a retry resumes rather than restarting.
BOLD_ATTEMPTS = 3
BOLD_RETRY_BACKOFF_SEC = 30


def _read_id_engine_hits(output_dir: Path, fasta_stem: str):
	"""Raw BOLD identification hits, from whichever store the installed version wrote.

	Both generations persist the hits — taxonomy, ``pct_identity``, ``process_id`` —
	before the fragile metadata phase runs, which is exactly what makes recovery
	possible when that phase dies:

	* 3.x — DuckDB at ``boldigger3_data/<stem>.duckdb``, table ``id_engine_results``.
	* 1.x — HDF at ``<stem>_result_storage.h5.lz``, key ``/results_unsorted``.
	"""
	if pd is None:
		return None

	for candidate in (
		output_dir / "boldigger3_data" / f"{fasta_stem}.duckdb",
		*sorted(output_dir.rglob("*.duckdb")),
	):
		if not candidate.exists() or candidate.stat().st_size == 0:
			continue
		try:
			import duckdb

			with duckdb.connect(str(candidate), read_only=True) as con:
				tables = [name for (name,) in con.execute("SHOW TABLES").fetchall()]
				# final_results is the metadata-joined table; id_engine_results is the
				# raw one that exists even when the metadata phase never ran.
				for table in ("final_results", "id_engine_results"):
					if table in tables:
						frame = con.execute(f"SELECT * FROM {table}").fetchdf()
						if frame is not None and not frame.empty:
							return frame
		except Exception:
			continue

	for candidate in (
		output_dir / f"{fasta_stem}_result_storage.h5.lz",
		*sorted(output_dir.rglob("*_result_storage.h5.lz")),
	):
		if not candidate.exists() or candidate.stat().st_size == 0:
			continue
		try:
			with pd.HDFStore(candidate, mode="r") as store:
				if "/results_unsorted" in store.keys():
					frame = store["/results_unsorted"]
					if frame is not None and not frame.empty:
						return frame
		except Exception:
			continue
	return None


# BOLD's own convention: how similar a hit must be before it is allowed to name a
# taxon at a given rank. 97% identity supports a species; 90% supports a family.
_BOLD_RANK_BY_THRESHOLD = ["Species", "Genus", "Family", "Order"]
_BOLD_RANK_ORDER = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def _match_rank_column(rank: str, columns) -> Optional[str]:
	"""The DataFrame column holding a rank, matched case-insensitively.

	BOLDigger3 3.x writes the Linnaean ranks LOWERCASE (``phylum``..``species``) in the
	``id_engine_results`` DuckDB, but ``_BOLD_RANK_ORDER`` is capitalised (for the output
	table and the case-insensitive consensus reader). Matching input columns by exact case
	found nothing and blanked every rank even for a 99% hit — a positive BOLD result would
	arrive with identity but no taxonomy. Match by folded case so the real lowercase store,
	and any capitalised variant, both resolve.
	"""
	target = rank.strip().lower()
	for column in columns:
		if str(column).strip().lower() == target:
			return str(column)
	return None


def select_bold_top_hits(hits, thresholds: List[int]) -> "pd.DataFrame":
	"""Collapse BOLD's many hits per query into one top hit each.

	Mirrors boldigger3's own ``select_top_hit``: the best identity in the group
	decides how deep the call is allowed to go, and the winning taxon at that rank is
	the most frequent one among the hits that clear the threshold. Ranks below the
	permitted level are left blank rather than guessed.
	"""
	rows = []
	for query_id, group in hits.groupby("id", sort=False):
		identities = group["pct_identity"].astype(float)
		top_identity = float(identities.max())

		level = "Class"
		cutoff = 0.0
		for index, rank in enumerate(_BOLD_RANK_BY_THRESHOLD):
			if index < len(thresholds) and top_identity >= float(thresholds[index]):
				level = rank
				cutoff = float(thresholds[index])
				break

		eligible = group[identities >= cutoff] if cutoff else group
		if eligible.empty:
			eligible = group

		row = {"id": query_id, "pct_identity": round(top_identity, 3)}
		allowed = _BOLD_RANK_ORDER[: _BOLD_RANK_ORDER.index(level) + 1]
		for rank in _BOLD_RANK_ORDER:
			source_col = _match_rank_column(rank, eligible.columns)
			if source_col is None or rank not in allowed:
				row[rank] = ""
				continue
			values = eligible[source_col].dropna().astype(str)
			values = values[values.str.strip().astype(bool)]
			row[rank] = values.mode().iloc[0] if not values.empty else ""
		for extra in ("bin_uri", "process_id"):
			if extra in eligible.columns and not eligible[extra].dropna().empty:
				row[extra] = str(eligible[extra].dropna().iloc[0])
		rows.append(row)
	return pd.DataFrame(rows)


def _salvage_bold_results(output_dir: Path, fasta_stem: str, config: BOLDConfig) -> dict[str, Path]:
	"""Rebuild BOLD's result tables from the identification hits already on disk.

	boldigger3 runs in three phases: identify (remote), join specimen metadata from
	the local DuckDB, then select the top hit. Both failure modes we hit — the
	``KeyError: 503`` when BOLD rate-limits, and the missing ``bold_public`` table —
	strike in phase two or three, *after* phase one has already persisted every hit
	with its full taxonomy and % identity.

	So the identification is never actually the thing that failed, and throwing it
	away would be the real loss. Re-deriving the top hit here recovers a complete,
	usable result; only cosmetic specimen metadata is missing.
	"""
	hits = _read_id_engine_hits(output_dir, fasta_stem)
	if hits is None or hits.empty:
		return {}
	if "id" not in hits.columns or "pct_identity" not in hits.columns:
		return {}
	hits = hits.dropna(subset=["id", "pct_identity"])
	if hits.empty:
		return {}

	files: dict[str, Path] = {}
	# Every hit, so the ranked hits table can offer the runner-up.
	all_hits_path = output_dir / "bold_all_hits.csv"
	hits.to_csv(all_hits_path, index=False)
	files["bold_all_hits_csv"] = all_hits_path

	top = select_bold_top_hits(hits, list(config.thresholds))
	if not top.empty:
		top_path = output_dir / "bold_top_hit.csv"
		top.to_csv(top_path, index=False)
		files["bold_parquet"] = top_path  # consumed by the consensus reader
	return files


def run_boldigger(fasta_path: Path, config: BOLDConfig, output_dir: Path) -> StepOutput:
	"""Run BOLDigger3 on a FASTA file and return a StepOutput."""
	output_dir.mkdir(parents=True, exist_ok=True)
	log_path = output_dir / "boldigger.log"
	sanitized_fasta = output_dir / "sanitized_input.fasta"
	header_mapping = _sanitize_fasta_for_bold(fasta_path, sanitized_fasta)

	errors = config.validate()
	if errors:
		log_path.write_text("\n".join(errors), encoding="utf-8")
		return StepOutput(step_name="boldigger", success=False, files={}, log=log_path, error_message="\n".join(errors))

	log_lines: list[str] = []

	# v3 needs a local DuckDB to exist and to hold a ``bold_public`` table, or its
	# metadata phase aborts and takes a successful identification down with it.
	version = boldigger_version()
	if version is not None and version >= MIN_BOLDIGGER_VERSION:
		_, note = ensure_bold_metadata_db(config.db_path)
		if note:
			log_lines.append(note)

	# BOLDigger3 drives the public BOLD v5 identification engine and does not log
	# in, so credentials are optional. They are still forwarded on stdin when
	# configured, for builds that do prompt.
	credentials = get_bold_credentials()
	credential_input = f"{credentials.username}\n{credentials.password}\n" if credentials else "\n"
	command = _boldigger_command(sanitized_fasta, config)

	last_error = ""
	for attempt in range(1, BOLD_ATTEMPTS + 1):
		log_lines.append(f"--- boldigger3 attempt {attempt}/{BOLD_ATTEMPTS}: {' '.join(command)}")
		try:
			completed = subprocess.run(
				command,
				input=credential_input,
				capture_output=True,
				text=True,
				timeout=config.timeout,
				cwd=str(output_dir),
				check=True,
			)
			log_lines.append((completed.stdout or "") + "\n" + (completed.stderr or ""))
			last_error = ""
			break
		except subprocess.TimeoutExpired:
			last_error = f"BOLDigger3 exceeded the timeout of {config.timeout}s."
			log_lines.append(last_error)
			break
		except subprocess.CalledProcessError as exc:
			output = (exc.stdout or "") + "\n" + (exc.stderr or "")
			log_lines.append(output)
			last_error = f"boldigger3 failed: {str(exc.stderr or exc.stdout or exc)[:400]}"
			# A 503 is BOLD rate-limiting us, not a bad request: waiting and resuming
			# is the correct response, and the run so far is preserved on disk.
			transient = "503" in output or "KeyError" in output
			if transient and attempt < BOLD_ATTEMPTS:
				log_lines.append(f"BOLD unavailable (503). Retrying in {BOLD_RETRY_BACKOFF_SEC}s...")
				time.sleep(BOLD_RETRY_BACKOFF_SEC)
				continue
			break

	files = _find_boldigger_outputs(output_dir, sanitized_fasta.stem)
	salvaged = False
	if not files:
		# boldigger3 never reached select_top_hit. Its identification results are
		# still on disk, so rebuild the top hits from them rather than losing the run.
		files = _salvage_bold_results(output_dir, sanitized_fasta.stem, config)
		if files:
			salvaged = True
			log_lines.append(
				"BOLDigger3 crashed after identifying (known bug: it does not handle BOLD's 503 code). "
				"Results recovered from the identification store; only the extra specimen metadata "
				"(sex, institution) was lost."
			)

	for path in files.values():
		_remap_bold_output_ids(path, header_mapping)
	log_path.write_text("\n".join(log_lines), encoding="utf-8")

	if not files:
		message = last_error or "BOLDigger3 exited cleanly but produced no outputs."
		return StepOutput(step_name="boldigger", success=False, files={}, log=log_path, error_message=message)

	return StepOutput(
		step_name="boldigger",
		success=True,
		files=files,
		log=log_path,
		error_message=("BOLD results salvaged after a boldigger3 failure (see log)." if salvaged else None),
	)


def download_bold_db(target_dir: Path) -> bool:
	"""Download the BOLD DuckDB as an explicit setup step.

	This is the step that needs a BOLD account: ``download_db`` asks for the
	username with ``input()``, the password with ``getpass()`` and then an extra
	``yes/no`` confirmation, whereas ``identify`` only reads the resulting local
	``*.ddb`` file.

	Best effort only. ``getpass()`` reads from ``/dev/tty`` when one is attached,
	so piping the password in cannot be relied on from an interactive shell —
	prefer running ``boldigger3 download_db <dir>`` by hand once.
	"""
	target_dir.mkdir(parents=True, exist_ok=True)
	credentials = get_bold_credentials()
	if credentials is None:
		return False
	stdin_text = f"{credentials.username}\n{credentials.password}\nyes\n"
	try:
		subprocess.run(
			[*_boldigger_base_command(), "download_db", str(target_dir)],
			input=stdin_text,
			capture_output=True,
			text=True,
			check=True,
			timeout=7200,
		)
		return find_local_bold_db(target_dir) is not None
	except (OSError, subprocess.SubprocessError):
		return False


def get_bold_reference(marker: str) -> str:
	"""Legacy symbolic helper preserved for compatibility."""
	_ = BOLD_CACHE_ROOT
	return f"bold:{marker.lower()}"