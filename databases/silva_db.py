"""SILVA reference database helpers."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from core.config import DEFAULT_SILVA_PATH, DEFAULT_SILVA_URL
from core.utils import (
	_gunzip_file,
	SINTAX_FORMAT_VERSION,
	load_silva_rank_map,
	silva_declared_labels,
	silva_rank_labels,
	taxonomy_fields_to_sintax,
)


def _convert_silva_to_sintax(src_fasta: Path, dest_fasta: Path, rank_map: Optional[dict] = None) -> None:
	"""Rewrite ARB-SILVA release headers into VSEARCH SINTAX format.

	SILVA headers look like ">accession.start.stop Domain;Phylum;Class;...";
	SINTAX needs ">accession;tax=d:Domain,p:Phylum,...;". Streams line by line
	since the SILVA release fasta can be several hundred MB.

	Without ``rank_map`` only the domain and the organism name are claimed
	(``silva_rank_labels``). With SILVA's rank file loaded into ``rank_map``, the
	intermediate levels get the rank SILVA declares for them (``silva_declared_labels``).
	If not a single header gains a rank from the map, the file is not the rank file of
	this release: the conversion is refused instead of quietly producing the unranked
	reference under a name that says it is ranked.
	"""
	dest_fasta.parent.mkdir(parents=True, exist_ok=True)
	tmp_dest = dest_fasta.with_suffix(dest_fasta.suffix + ".tmp")
	headers_with_declared_rank = 0
	with src_fasta.open("r", encoding="utf-8", errors="replace") as src, tmp_dest.open("w", encoding="utf-8") as dst:
		for line in src:
			if line.startswith(">"):
				header = line[1:].rstrip("\n")
				accession, _, taxonomy = header.partition(" ")
				fields = [f for f in taxonomy.split(";") if f.strip()]
				# SILVA lineages vary in depth and always end with the organism name,
				# so the ranks must be anchored rather than filled positionally.
				if rank_map:
					labels = silva_declared_labels(fields, rank_map)
					if any(label not in (None, "d", "s") for label in labels):
						headers_with_declared_rank += 1
				else:
					labels = silva_rank_labels(len(fields))
				dst.write(f">{accession};tax={taxonomy_fields_to_sintax(fields, labels)};\n")
			else:
				dst.write(line)
	if rank_map is not None and headers_with_declared_rank == 0:
		tmp_dest.unlink(missing_ok=True)
		raise ValueError(
			"the SILVA rank file matched no lineage in the SILVA FASTA "
			"(is it the tax_slv_ssu file of the same release?)"
		)
	tmp_dest.replace(dest_fasta)


def _file_digest(path: Path, length: int = 8) -> str:
	"""Short sha256 of a file's content, to name a cache after the exact input it used."""
	digest = hashlib.sha256()
	with Path(path).open("rb") as handle:
		for chunk in iter(lambda: handle.read(1 << 20), b""):
			digest.update(chunk)
	return digest.hexdigest()[:length]


def silva_rank_path(db_paths) -> tuple:
	"""``(path, note)``: the SILVA rank file to use, and a log line when it cannot be.

	``path`` is ``None`` when ``DB_SILVA_RANKS_PATH`` is unset (no note: the default) or
	set to a file that does not exist (with a note saying so, so a typo in .env shows up
	in the run log instead of silently giving the unranked reference).
	"""
	configured = getattr(db_paths, "silva_ranks", None)
	if not configured:
		return None, None
	configured = Path(configured)
	if not configured.exists():
		return None, (
			f"DB_SILVA_RANKS_PATH points to a missing file ({configured}); "
			"SILVA runs without declared ranks (domain and organism name only)."
		)
	return configured, None


@dataclass(frozen=True)
class SILVAConfig:
	db_path: Path
	version: str = "138.2"
	target: str = "SSU"
	# SILVA's rank file. None keeps the original behaviour (domain + organism name).
	rank_path: Optional[Path] = None

	def validate(self) -> list[str]:
		errors = []
		if not self.db_path.exists():
			errors.append(f"SILVA database not found: {self.db_path}")
		valid_suffixes = {".fasta", ".fa", ".fna", ".gz"}
		if self.db_path.suffix not in valid_suffixes:
			errors.append(f"Unrecognised SILVA file extension: {self.db_path.suffix}")
		return errors

	def _decompressed_path(self) -> Path:
		if self.db_path.suffix.lower() != ".gz":
			return self.db_path
		decompressed = self.db_path.with_suffix("")
		if not decompressed.exists() or decompressed.stat().st_mtime < self.db_path.stat().st_mtime:
			_gunzip_file(str(self.db_path), str(decompressed))
		return decompressed

	def sintax_ready_path(self) -> Path:
		"""Return a SINTAX-formatted copy of this SILVA fasta, converting/caching on demand."""
		working = self._decompressed_path()
		if self.rank_path is None:
			# Versioned so a cache built by the old positional rank mapping is rebuilt.
			# Unchanged name: references built before the rank file existed are reused.
			sintax_path = working.parent / f"{working.stem}.sintax.{SINTAX_FORMAT_VERSION}.fasta"
			if not sintax_path.exists() or sintax_path.stat().st_mtime < working.stat().st_mtime:
				_convert_silva_to_sintax(working, sintax_path)
			return sintax_path
		# The ranked copy lives under its own name, keyed by the rank file's content, so
		# it never overwrites the unranked one (runs made with it stay reproducible) and a
		# different rank file can never be served a cache built from another.
		rank_tag = _file_digest(self.rank_path)
		sintax_path = working.parent / f"{working.stem}.sintax.{SINTAX_FORMAT_VERSION}r.{rank_tag}.fasta"
		if not sintax_path.exists() or sintax_path.stat().st_mtime < working.stat().st_mtime:
			rank_map = load_silva_rank_map(self.rank_path)
			if not rank_map:
				raise ValueError(f"not a SILVA rank file (no path/rank lines): {self.rank_path}")
			_convert_silva_to_sintax(working, sintax_path, rank_map=rank_map)
		return sintax_path

	def rank_note(self) -> str:
		"""One line for logs and source_status: which rank file the reference used."""
		if self.rank_path is None:
			return ""
		return f"declared ranks from {Path(self.rank_path).name}"

	def as_sintax_args(self) -> list[str]:
		return ["--db", str(self.sintax_ready_path())]

	def as_dada2_ref(self) -> str:
		return str(self.db_path)


def load_silva(db_path: Path, version: str = "138.2", target: str = "SSU") -> SILVAConfig:
	config = SILVAConfig(db_path=db_path, version=version, target=target)
	errors = config.validate()
	if errors:
		raise ValueError(f"SILVA invalid: {'; '.join(errors)}")
	return config


def get_silva_reference(marker: str) -> str:
	"""Legacy symbolic helper preserved for compatibility."""
	_ = (DEFAULT_SILVA_PATH, DEFAULT_SILVA_URL)
	return f"silva:{marker.lower()}"