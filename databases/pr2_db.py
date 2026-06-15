"""PR2 reference database helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from core.config import DEFAULT_PR2_DADA2_FASTA_PATH, DEFAULT_PR2_MOTHUR_FASTA_PATH, DEFAULT_PR2_MOTHUR_TAX_PATH, DEFAULT_PR2_URL


@dataclass(frozen=True)
class PR2Config:
	db_path: Path
	version: str = "5.0"

	def validate(self) -> list[str]:
		errors = []
		if not self.db_path.exists():
			errors.append(f"Base PR2 nao encontrada: {self.db_path}")
		if str(self.db_path).lower().endswith(".gz"):
			errors.append(
				f"PR2 comprimida nao pode ser usada directamente com VSEARCH SINTAX: {self.db_path}. "
				"Descomprime primeiro."
			)
		return errors

	def as_sintax_args(self) -> list[str]:
		return ["--db", str(self.db_path)]

	def as_dada2_ref(self) -> str:
		return str(self.db_path)


def load_pr2(db_path: Path, version: str = "5.0") -> PR2Config:
	config = PR2Config(db_path=db_path, version=version)
	errors = config.validate()
	if errors:
		raise ValueError(f"PR2 invalido: {'; '.join(errors)}")
	return config


def _normalize_pr2_mothur_db(fasta_in, tax_in, out_fasta, out_map):
	_ = (fasta_in, tax_in, out_fasta, out_map)
	return 0


def prepare_pr2_mothur_db(pr2_dir):
	_ = (DEFAULT_PR2_URL, DEFAULT_PR2_MOTHUR_FASTA_PATH, DEFAULT_PR2_MOTHUR_TAX_PATH, DEFAULT_PR2_DADA2_FASTA_PATH)
	return pr2_dir


def get_pr2_reference(marker: str) -> str:
	"""Legacy symbolic helper preserved for compatibility."""
	return f"pr2:{marker.lower()}"