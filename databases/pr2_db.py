"""PR2 reference database helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from core.config import DEFAULT_PR2_DADA2_FASTA_PATH, DEFAULT_PR2_MOTHUR_FASTA_PATH, DEFAULT_PR2_MOTHUR_TAX_PATH, DEFAULT_PR2_URL
from core.utils import _gunzip_file, is_pr2_placeholder, PR2_RANK_LABELS, SINTAX_FORMAT_VERSION, taxonomy_fields_to_sintax


def _normalize_pr2_mothur_db(fasta_in: Path, tax_in: Path, out_fasta: Path) -> int:
	"""Join PR2's mothur-flavor fasta (headers are bare IDs) with its separate
	tab-separated taxonomy file into a single SINTAX-ready fasta. Returns the
	number of fasta records that had no matching taxonomy entry (skipped).
	"""
	tax_by_id = {}
	with tax_in.open("r", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			line = line.rstrip("\n")
			if not line:
				continue
			seq_id, _, taxonomy = line.partition("\t")
			# Blank out PR2's unresolved-lineage placeholders (Craniata_X, *_sp.) rather
			# than drop them: the nine levels are positional, and removing one would
			# shift every level below it onto the wrong rank. An empty field is skipped
			# by taxonomy_fields_to_sintax without disturbing its neighbours.
			fields = [
				"" if is_pr2_placeholder(f) else f.strip()
				for f in taxonomy.strip().strip(";").split(";")
			]
			tax_by_id[seq_id] = fields

	out_fasta.parent.mkdir(parents=True, exist_ok=True)
	tmp_out = out_fasta.with_suffix(out_fasta.suffix + ".tmp")
	skipped = 0
	with fasta_in.open("r", encoding="utf-8", errors="replace") as src, tmp_out.open("w", encoding="utf-8") as dst:
		write_seq = False
		for line in src:
			if line.startswith(">"):
				seq_id = line[1:].strip().split()[0]
				fields = tax_by_id.get(seq_id)
				if fields is None:
					skipped += 1
					write_seq = False
					continue
				# PR2's nine levels are not Linnean; without its explicit rank map a
				# family lands on the species rank (see PR2_RANK_LABELS).
				dst.write(f">{seq_id};tax={taxonomy_fields_to_sintax(fields, PR2_RANK_LABELS)};\n")
				write_seq = True
			elif write_seq:
				dst.write(line)
	tmp_out.replace(out_fasta)
	return skipped


@dataclass(frozen=True)
class PR2Config:
	db_path: Path
	tax_path: Optional[Path] = None
	version: str = "5.0"

	def validate(self) -> list[str]:
		errors = []
		if not self.db_path.exists():
			errors.append(f"PR2 database not found: {self.db_path}")
		return errors

	def _decompressed(self, path: Path) -> Path:
		if path.suffix.lower() != ".gz":
			return path
		decompressed = path.with_suffix("")
		if not decompressed.exists() or decompressed.stat().st_mtime < path.stat().st_mtime:
			_gunzip_file(str(path), str(decompressed))
		return decompressed

	def sintax_ready_path(self) -> Path:
		"""Return a SINTAX-formatted copy of this PR2 mothur fasta+tax pair, converting/caching on demand."""
		if not self.tax_path or not self.tax_path.exists():
			raise FileNotFoundError(f"PR2 taxonomy file (tax_path) not configured or not found: {self.tax_path}")
		fasta_working = self._decompressed(self.db_path)
		tax_working = self._decompressed(self.tax_path)
		# The cache name carries the conversion's format version: mtime alone cannot
		# tell that the *rank mapping* changed, so a stale cache from before
		# PR2_RANK_LABELS would keep feeding families in as species.
		sintax_path = fasta_working.parent / f"{fasta_working.stem}.sintax.{SINTAX_FORMAT_VERSION}.fasta"
		if (
			not sintax_path.exists()
			or sintax_path.stat().st_mtime < fasta_working.stat().st_mtime
			or sintax_path.stat().st_mtime < tax_working.stat().st_mtime
		):
			_normalize_pr2_mothur_db(fasta_working, tax_working, sintax_path)
		return sintax_path

	def as_sintax_args(self) -> list[str]:
		return ["--db", str(self.sintax_ready_path())]

	def as_dada2_ref(self) -> str:
		return str(self.db_path)


def load_pr2(db_path: Path, version: str = "5.0") -> PR2Config:
	config = PR2Config(db_path=db_path, version=version)
	errors = config.validate()
	if errors:
		raise ValueError(f"PR2 invalid: {'; '.join(errors)}")
	return config


def prepare_pr2_mothur_db(pr2_dir):
	_ = (DEFAULT_PR2_URL, DEFAULT_PR2_MOTHUR_FASTA_PATH, DEFAULT_PR2_MOTHUR_TAX_PATH, DEFAULT_PR2_DADA2_FASTA_PATH)
	return pr2_dir


def get_pr2_reference(marker: str) -> str:
	"""Legacy symbolic helper preserved for compatibility."""
	return f"pr2:{marker.lower()}"