"""SILVA reference database helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from core.config import DEFAULT_SILVA_PATH, DEFAULT_SILVA_URL


@dataclass(frozen=True)
class SILVAConfig:
	db_path: Path
	version: str = "138.2"
	target: str = "SSU"

	def validate(self) -> list[str]:
		errors = []
		if not self.db_path.exists():
			errors.append(f"Base SILVA nao encontrada: {self.db_path}")
		if str(self.db_path).lower().endswith(".gz"):
			errors.append(
				f"SILVA comprimida nao pode ser usada directamente com VSEARCH SINTAX: {self.db_path}. "
				"Descomprime primeiro com: gunzip silva.fasta.gz"
			)
		valid_suffixes = {".fasta", ".fa", ".fna", ".gz"}
		if self.db_path.suffix not in valid_suffixes:
			errors.append(f"Extensao nao reconhecida para SILVA: {self.db_path.suffix}")
		return errors

	def as_sintax_args(self) -> list[str]:
		return ["--db", str(self.db_path)]

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