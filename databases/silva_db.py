"""SILVA reference database helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from core.config import DEFAULT_SILVA_PATH, DEFAULT_SILVA_URL
from core.utils import _gunzip_file, SINTAX_FORMAT_VERSION, silva_rank_labels, taxonomy_fields_to_sintax


def _convert_silva_to_sintax(src_fasta: Path, dest_fasta: Path) -> None:
	"""Rewrite ARB-SILVA release headers into VSEARCH SINTAX format.

	SILVA headers look like ">accession.start.stop Domain;Phylum;Class;...";
	SINTAX needs ">accession;tax=d:Domain,p:Phylum,...;". Streams line by line
	since the SILVA release fasta can be several hundred MB.
	"""
	dest_fasta.parent.mkdir(parents=True, exist_ok=True)
	tmp_dest = dest_fasta.with_suffix(dest_fasta.suffix + ".tmp")
	with src_fasta.open("r", encoding="utf-8", errors="replace") as src, tmp_dest.open("w", encoding="utf-8") as dst:
		for line in src:
			if line.startswith(">"):
				header = line[1:].rstrip("\n")
				accession, _, taxonomy = header.partition(" ")
				fields = [f for f in taxonomy.split(";") if f.strip()]
				# SILVA lineages vary in depth and always end with the organism name,
				# so the ranks must be anchored rather than filled positionally.
				dst.write(f">{accession};tax={taxonomy_fields_to_sintax(fields, silva_rank_labels(len(fields)))};\n")
			else:
				dst.write(line)
	tmp_dest.replace(dest_fasta)


@dataclass(frozen=True)
class SILVAConfig:
	db_path: Path
	version: str = "138.2"
	target: str = "SSU"

	def validate(self) -> list[str]:
		errors = []
		if not self.db_path.exists():
			errors.append(f"Base SILVA nao encontrada: {self.db_path}")
		valid_suffixes = {".fasta", ".fa", ".fna", ".gz"}
		if self.db_path.suffix not in valid_suffixes:
			errors.append(f"Extensao nao reconhecida para SILVA: {self.db_path.suffix}")
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
		# Versioned so a cache built by the old positional rank mapping is rebuilt.
		sintax_path = working.parent / f"{working.stem}.sintax.{SINTAX_FORMAT_VERSION}.fasta"
		if not sintax_path.exists() or sintax_path.stat().st_mtime < working.stat().st_mtime:
			_convert_silva_to_sintax(working, sintax_path)
		return sintax_path

	def as_sintax_args(self) -> list[str]:
		return ["--db", str(self.sintax_ready_path())]

	def as_dada2_ref(self) -> str:
		return str(self.db_path)


def load_silva(db_path: Path, version: str = "138.2", target: str = "SSU") -> SILVAConfig:
	config = SILVAConfig(db_path=db_path, version=version, target=target)
	errors = config.validate()
	if errors:
		raise ValueError(f"SILVA invalido: {'; '.join(errors)}")
	return config


def get_silva_reference(marker: str) -> str:
	"""Legacy symbolic helper preserved for compatibility."""
	_ = (DEFAULT_SILVA_PATH, DEFAULT_SILVA_URL)
	return f"silva:{marker.lower()}"