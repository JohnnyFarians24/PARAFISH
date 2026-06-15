"""Taxonomy classification helpers for BOLDigger, SINTAX and consensus building."""

from __future__ import annotations

import csv
import subprocess
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

try:
	import pandas as pd
except Exception:  # pragma: no cover
	pd = None

from core.config import BOLDIGGER_PYTHON, BOLDIGGER_TIMEOUT_SEC, VSEARCH_PATH
from core.models import StepOutput

from databases.bold_db import BOLDConfig, run_boldigger as run_boldigger_step
from databases.ncbi_db import NCBIConfig, run_blast
from databases.pr2_db import PR2Config
from databases.silva_db import SILVAConfig

def _run_boldigger_legacy(
	fasta_path: Path,
	bold_params: Dict,
	db_path: Path,
	output_dir: Optional[Path] = None,
	prefix: str = "bold",
) -> Path:
	"""Compatibility shim that delegates to the real BOLDigger runner.

	This preserves the historical return type used by older callers while
	avoiding the old behavior of overwriting real BOLDigger outputs.
	"""
	warnings.warn(
		"databases.taxonomy_runner.run_boldigger is deprecated; use databases.bold_db.run_boldigger_step instead.",
		DeprecationWarning,
		stacklevel=2,
	)
	output_dir = output_dir or (fasta_path.parent / "taxonomy")
	output_dir.mkdir(parents=True, exist_ok=True)
	config = BOLDConfig(
		db_path=Path(db_path),
		db_source=int((bold_params or {}).get("db", 1)),
		mode=int((bold_params or {}).get("mode", 2)),
		thresholds=[int(value) for value in (bold_params or {}).get("thresholds", [97, 95, 90, 85])],
	)
	result = run_boldigger_step(fasta_path, config, output_dir)
	if result.files:
		return next(iter(result.files.values()))
	fallback = output_dir / f"{prefix}_result.parquet"
	if not fallback.exists():
		fallback.touch()
	return fallback


run_boldigger = _run_boldigger_legacy


TAXONOMY_LEVELS = ["kingdom", "phylum", "class_", "order", "family", "genus", "species"]


class ConsensusTable:
	"""Lightweight table container used when pandas is unavailable."""

	def __init__(self, rows: List[Dict[str, str]], columns: List[str]):
		self.rows = rows
		self.columns = columns

	def to_csv(self, path: Path) -> None:
		path.parent.mkdir(parents=True, exist_ok=True)
		with path.open("w", newline="", encoding="utf-8") as handle:
			writer = csv.DictWriter(handle, fieldnames=self.columns)
			writer.writeheader()
			writer.writerows(self.rows)

	def __len__(self):
		return len(self.rows)


def _read_table_rows(path: Path) -> List[Dict[str, str]]:
	if not path or not path.exists() or path.stat().st_size == 0:
		return []
	if path.suffix.lower() in {".parquet", ".xlsx"}:
		try:
			if pd is None:
				raise ImportError("pandas not available")
			if path.suffix.lower() == ".parquet":
				df = pd.read_parquet(path)
			else:
				df = pd.read_excel(path, engine="openpyxl")
			df = df.where(pd.notnull(df), None)
			return [dict(row) for row in df.to_dict(orient="records")]
		except Exception:
			with path.open("r", encoding="utf-8", errors="replace") as handle:
				return list(csv.DictReader(handle))
	with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
		sample = handle.read(4096)
		handle.seek(0)
		dialect = csv.Sniffer().sniff(sample) if sample.strip() else csv.excel
		return list(csv.DictReader(handle, dialect=dialect))


def _safe_text(value) -> str:
	return str(value).strip() if value is not None else ""


def _row_value(row: Dict[str, str], key: str) -> str:
	for candidate in (key, key.lower(), key.upper()):
		if candidate in row:
			return _safe_text(row.get(candidate))
	return ""


def _normalize_taxon_value(value: str) -> str:
	clean = value.strip()
	if not clean:
		return ""
	if clean.lower() in {"na", "unclassified", "unknown", "none", "nan"}:
		return ""
	return clean


def _first_nonempty(row: Dict[str, str], *keys: str) -> str:
	for key in keys:
		value = _row_value(row, key)
		if value:
			return value
	return ""


def _extract_asv_id(row: Dict[str, str]) -> str:
	raw = _first_nonempty(row, "ASV_ID", "asv_id", "qseqid", "id", "seq_id", "query")
	return raw.split(";")[0].strip() if raw else ""


def _normalize_source_row(source: str, row: Dict[str, str]) -> Dict[str, str]:
	normalized = {
		"asv_id": _extract_asv_id(row),
		"source": source,
		"kingdom": _first_nonempty(row, "kingdom", "Consensus_Kingdom"),
		"phylum": _first_nonempty(row, "phylum", "Consensus_Phylum"),
		"class_": _first_nonempty(row, "class_", "class", "Consensus_Class"),
		"order": _first_nonempty(row, "order", "Consensus_Order"),
		"family": _first_nonempty(row, "family", "Consensus_Family"),
		"genus": _first_nonempty(row, "genus", "Consensus_Genus"),
		"species": _first_nonempty(row, "species", "Consensus_Species", "SINTAX_species", "BOLD_species", "sci_name", "sscinames"),
		"similarity": _first_nonempty(row, "similarity", "BOLD_similarity", "pident"),
		"confidence": _first_nonempty(row, "confidence", "SINTAX_bootstrap", "bootstrap", "bitscore"),
		"flags": _first_nonempty(row, "flags", "BOLD_flags", "Warnings"),
	}
	return normalized


def _lca_consensus(records: list[dict]) -> dict:
	"""Return the deepest shared taxonomy level across records."""
	consensus: dict = {}
	if not records:
		return consensus
	for level in reversed(TAXONOMY_LEVELS):
		values = [
			_normalize_taxon_value(_safe_text(record.get(level))).lower()
			for record in records
			if _normalize_taxon_value(_safe_text(record.get(level)))
		]
		values = [value for value in values if value]
		if not values:
			continue
		if len(set(values)) == 1:
			original = next(
				_safe_text(record.get(level))
				for record in records
				if _normalize_taxon_value(_safe_text(record.get(level)))
			)
			consensus[level] = original
			idx = TAXONOMY_LEVELS.index(level)
			for upper_level in TAXONOMY_LEVELS[:idx]:
				candidates = [_safe_text(record.get(upper_level)) for record in records if _safe_text(record.get(upper_level))]
				if candidates:
					consensus[upper_level] = candidates[0]
			break
	if not consensus:
		kingdom_values = {_normalize_taxon_value(_safe_text(record.get("kingdom"))) for record in records if _safe_text(record.get("kingdom"))}
		kingdom_values = {value for value in kingdom_values if value}
		if len(kingdom_values) == 1:
			consensus["kingdom"] = next(iter(kingdom_values))
		else:
			consensus["kingdom"] = "unclassified"
	return consensus


def _strict_consensus(records: list[dict]) -> dict:
	"""Return the deepest taxonomy path identical across every record."""
	consensus: dict = {}
	if not records:
		return consensus
	for level in TAXONOMY_LEVELS:
		values = [_normalize_taxon_value(_safe_text(record.get(level))) for record in records]
		values = [value for value in values if value]
		if not values:
			continue
		if len({value.lower() for value in values}) == 1:
			consensus[level] = values[0]
			continue
		break
	if not consensus:
		return _lca_consensus(records)
	return consensus


def _majority_consensus(records: list[dict]) -> dict:
	"""Return the deepest taxonomy level supported by a strict majority."""
	consensus: dict = {}
	if not records:
		return consensus
	for level in reversed(TAXONOMY_LEVELS):
		values = [_normalize_taxon_value(_safe_text(record.get(level))) for record in records]
		values = [value for value in values if value]
		if not values:
			continue
		counts: Dict[str, int] = defaultdict(int)
		originals: Dict[str, str] = {}
		for value in values:
			key = value.lower()
			counts[key] += 1
			originals.setdefault(key, value)
		best_key, best_count = max(counts.items(), key=lambda item: (item[1], item[0]))
		if best_count > len(records) / 2:
			consensus[level] = originals[best_key]
			idx = TAXONOMY_LEVELS.index(level)
			for upper_level in TAXONOMY_LEVELS[:idx]:
				candidates = [_safe_text(record.get(upper_level)) for record in records if _safe_text(record.get(upper_level))]
				if candidates:
					consensus[upper_level] = candidates[0]
			break
	if not consensus:
		return _lca_consensus(records)
	return consensus


def _build_consensus_row(asv_id: str, records: list[dict], method: str) -> Dict[str, str]:
	method_name = (method or "lca").strip().lower()
	if method_name == "strict":
		consensus = _strict_consensus(records)
	elif method_name == "majority":
		consensus = _majority_consensus(records)
	else:
		consensus = _lca_consensus(records)
	row = {
		"ASV_ID": asv_id,
		"Consensus_Kingdom": consensus.get("kingdom", ""),
		"Consensus_Phylum": consensus.get("phylum", ""),
		"Consensus_Class": consensus.get("class_", ""),
		"Consensus_Order": consensus.get("order", ""),
		"Consensus_Family": consensus.get("family", ""),
		"Consensus_Genus": consensus.get("genus", ""),
		"Consensus_Species": consensus.get("species", ""),
		"BOLD_species": "",
		"BOLD_similarity": "",
		"BOLD_flags": "",
		"SINTAX_species": "",
		"SINTAX_bootstrap": "",
		"Source_of_consensus": method_name,
		"Warnings": "",
	}
	for record in records:
		if record.get("source") == "BOLD":
			row["BOLD_species"] = record.get("species", "")
			row["BOLD_similarity"] = record.get("similarity", "")
			row["BOLD_flags"] = record.get("flags", "")
		elif record.get("source") in {"SILVA", "PR2", "CUSTOM"}:
			row["SINTAX_species"] = record.get("species", "")
			row["SINTAX_bootstrap"] = record.get("confidence", "")
	return row


def _annotate_fasta(fasta_path: Path, consensus_rows: Dict[str, Dict[str, str]], output_path: Path) -> None:
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with fasta_path.open("r", encoding="utf-8", errors="replace") as input_handle, output_path.open("w", encoding="utf-8") as output_handle:
		current_header = None
		sequence_lines: list[str] = []
		for raw_line in input_handle:
			line = raw_line.rstrip("\n")
			if line.startswith(">"):
				if current_header is not None:
					asv_id = current_header.split()[0]
					taxonomy = consensus_rows.get(asv_id, {})
					size = 1
					if "size=" in current_header:
						try:
							size = int(current_header.split("size=", 1)[1].split(";", 1)[0])
						except Exception:
							size = 1
					tax_parts = [taxonomy.get(key, "") for key in ["Consensus_Kingdom", "Consensus_Phylum", "Consensus_Class", "Consensus_Order", "Consensus_Family", "Consensus_Genus", "Consensus_Species"]]
					output_handle.write(f">{asv_id};size={size};tax={';'.join(tax_parts)}\n")
					for seq_line in sequence_lines:
						output_handle.write(seq_line + "\n")
				current_header = line[1:].strip()
				sequence_lines = []
			else:
				sequence_lines.append(line)
		if current_header is not None:
			asv_id = current_header.split()[0]
			taxonomy = consensus_rows.get(asv_id, {})
			size = 1
			if "size=" in current_header:
				try:
					size = int(current_header.split("size=", 1)[1].split(";", 1)[0])
				except Exception:
					size = 1
			tax_parts = [taxonomy.get(key, "") for key in ["Consensus_Kingdom", "Consensus_Phylum", "Consensus_Class", "Consensus_Order", "Consensus_Family", "Consensus_Genus", "Consensus_Species"]]
			output_handle.write(f">{asv_id};size={size};tax={';'.join(tax_parts)}\n")
			for seq_line in sequence_lines:
				output_handle.write(seq_line + "\n")


def _build_conflict_report(records_by_asv: Dict[str, List[dict]]) -> List[Dict[str, str]]:
	report_rows: List[Dict[str, str]] = []
	for asv_id, records in records_by_asv.items():
		if len(records) < 2:
			continue
		conflict_levels = []
		for level in ["genus", "family", "order", "class_", "phylum", "kingdom"]:
			values = {_normalize_taxon_value(_safe_text(record.get(level))).lower() for record in records if _safe_text(record.get(level))}
			values = {value for value in values if value}
			if len(values) > 1:
				conflict_levels.append(level)
		if conflict_levels:
			report_rows.append(
				{
					"ASV_ID": asv_id,
					"conflict_levels": ",".join(conflict_levels),
					"sources": ",".join(sorted({record.get("source", "") for record in records if record.get("source")})),
				}
			)
	return report_rows


def _write_csv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
	path.parent.mkdir(parents=True, exist_ok=True)
	with path.open("w", newline="", encoding="utf-8") as handle:
		writer = csv.DictWriter(handle, fieldnames=columns)
		writer.writeheader()
		writer.writerows(rows)


def run_sintax(fasta_path: Path, db_path: Path, cutoff: float, threads: int, output_dir: Optional[Path] = None, prefix: str = "sintax") -> Path:
	"""Run VSEARCH SINTAX and return the result CSV path."""
	output_dir = output_dir or (fasta_path.parent / "taxonomy")
	output_dir.mkdir(parents=True, exist_ok=True)
	result_path = output_dir / f"{prefix}_result.csv"
	log_path = output_dir / f"{prefix}.log"

	command = [
		VSEARCH_PATH,
		"--sintax",
		str(fasta_path),
		"--db",
		str(db_path),
		"--sintax_cutoff",
		str(cutoff),
		"--tabbedout",
		str(result_path),
		"--threads",
		str(threads),
		"--log",
		str(log_path),
	]
	try:
		completed = subprocess.run(command, capture_output=True, text=True, check=True)
		log_path.write_text((completed.stdout or "") + ("\n" if completed.stdout and completed.stderr else "") + (completed.stderr or ""), encoding="utf-8")
	except subprocess.CalledProcessError as exc:
		log_path.write_text(str(exc), encoding="utf-8")
	if not result_path.exists():
		result_path.write_text("ASV_ID,Taxonomy\n", encoding="utf-8")
	return result_path


def build_consensus(
	bold_result: Optional[Path],
	sintax_result: Optional[Path],
	blast_result: Optional[Path],
	method: str = "lca",
) -> ConsensusTable:
	"""Normalize and combine taxonomy results from multiple sources."""
	rows_by_asv: Dict[str, List[dict]] = defaultdict(list)
	for source, table_path in (("BOLD", bold_result), ("SILVA", sintax_result), ("NCBI", blast_result)):
		if table_path and table_path.exists():
			for row in _read_table_rows(table_path):
				normalized = _normalize_source_row(source, row)
				if normalized["asv_id"]:
					rows_by_asv[normalized["asv_id"]].append(normalized)
	consensus_rows = [_build_consensus_row(asv_id, records, method) for asv_id, records in sorted(rows_by_asv.items())]
	columns = [
		"ASV_ID",
		"Consensus_Kingdom",
		"Consensus_Phylum",
		"Consensus_Class",
		"Consensus_Order",
		"Consensus_Family",
		"Consensus_Genus",
		"Consensus_Species",
		"BOLD_species",
		"BOLD_similarity",
		"BOLD_flags",
		"SINTAX_species",
		"SINTAX_bootstrap",
		"Source_of_consensus",
		"Warnings",
	]
	return ConsensusTable(consensus_rows, columns)


def _collect_source_records(source: str, table_path: Optional[Path]) -> Dict[str, List[dict]]:
	rows_by_asv: Dict[str, List[dict]] = defaultdict(list)
	if table_path and table_path.exists() and table_path.stat().st_size > 0:
		for row in _read_table_rows(table_path):
			normalized = _normalize_source_row(source, row)
			if normalized["asv_id"]:
				rows_by_asv[normalized["asv_id"]].append(normalized)
	return rows_by_asv


def run_taxonomy_pipeline(
	fasta_path: Path,
	output_dir: Path,
	silva_config: Optional[SILVAConfig] = None,
	pr2_config: Optional[PR2Config] = None,
	ncbi_config: Optional[NCBIConfig] = None,
	bold_config: Optional[BOLDConfig] = None,
	custom_db_path: Optional[Path] = None,
	consensus_method: str = "lca",
	threads: int = 8,
	sintax_cutoff: float = 0.8,
) -> StepOutput:
	"""Run all active taxonomy sources and return a single consensus StepOutput."""
	output_dir.mkdir(parents=True, exist_ok=True)
	log_path = output_dir / "taxonomy_pipeline.log"
	consensus_path = output_dir / "taxonomy_consensus.csv"
	annotated_fasta = output_dir / "annotated.fasta"
	conflict_report = output_dir / "conflict_report.csv"
	files: Dict[str, Path] = {"log": log_path}
	metrics: Dict[str, float] = {}
	log_lines: List[str] = []
	records_by_asv: Dict[str, List[dict]] = defaultdict(list)
	active_sources = 0

	def _merge_records(source_name: str, table_path: Optional[Path]) -> None:
		nonlocal active_sources
		rows = _collect_source_records(source_name, table_path)
		if rows:
			active_sources += 1
		for asv_id, source_rows in rows.items():
			records_by_asv[asv_id].extend(source_rows)

	if silva_config is not None:
		silva_path = run_sintax(fasta_path, silva_config.db_path, sintax_cutoff, threads, output_dir=output_dir / "silva", prefix="silva_sintax")
		files["silva_sintax_csv"] = silva_path
		_merge_records("SILVA", silva_path)

	if pr2_config is not None:
		pr2_path = run_sintax(fasta_path, pr2_config.db_path, sintax_cutoff, threads, output_dir=output_dir / "pr2", prefix="pr2_sintax")
		files["pr2_sintax_csv"] = pr2_path
		_merge_records("PR2", pr2_path)

	if custom_db_path is not None and custom_db_path.exists():
		custom_path = run_sintax(fasta_path, custom_db_path, sintax_cutoff, threads, output_dir=output_dir / "custom", prefix="custom_sintax")
		files["custom_sintax_csv"] = custom_path
		_merge_records("CUSTOM", custom_path)

	if ncbi_config is not None:
		ncbi_result = run_blast(fasta_path, ncbi_config, output_dir / "ncbi")
		if ncbi_result.success:
			files.update(ncbi_result.files)
			blast_tsv = ncbi_result.files.get("blast_tsv")
			if isinstance(blast_tsv, Path) and blast_tsv.exists():
				_merge_records("NCBI", blast_tsv)
			else:
				metrics["ncbi_error"] = 1.0
				log_lines.append("NCBI result missing blast_tsv path.")
		else:
			metrics["ncbi_error"] = 1.0
			log_lines.append(ncbi_result.error_message or "NCBI failed")

	if bold_config is not None:
		ddb_files = list(bold_config.db_path.glob("*.ddb")) if bold_config.db_path.exists() else []
		if not ddb_files:
			log_lines.append(f"BOLD database empty or missing at {bold_config.db_path} — skipping BOLD.")
			bold_config = None

	if bold_config is not None:
		bold_result = run_boldigger_step(fasta_path, bold_config, output_dir / "bold")
		if bold_result.success:
			files.update(bold_result.files)
			bold_path = bold_result.files.get("bold_parquet") or bold_result.files.get("bold_xlsx")
			if isinstance(bold_path, Path) and bold_path.exists():
				_merge_records("BOLD", bold_path)
			else:
				metrics["bold_error"] = 1.0
				log_lines.append("BOLD result missing parquet/xlsx path.")
		else:
			metrics["bold_error"] = 1.0
			log_lines.append(bold_result.error_message or "BOLD failed")

	if not records_by_asv:
		log_path.write_text("\n".join(log_lines) or "No taxonomy sources succeeded.", encoding="utf-8")
		return StepOutput(step_name="taxonomy", success=False, files=files, log=log_path, metrics=metrics, error_message="Nenhuma fonte de taxonomia teve sucesso.")

	consensus_rows = [_build_consensus_row(asv_id, records, consensus_method) for asv_id, records in sorted(records_by_asv.items())]
	conflicts = _build_conflict_report(records_by_asv)
	conflict_columns = ["ASV_ID", "conflict_levels", "sources"]
	_write_csv(conflict_report, conflicts, conflict_columns)
	consensus_table = ConsensusTable(
		consensus_rows,
		[
			"ASV_ID",
			"Consensus_Kingdom",
			"Consensus_Phylum",
			"Consensus_Class",
			"Consensus_Order",
			"Consensus_Family",
			"Consensus_Genus",
			"Consensus_Species",
			"BOLD_species",
			"BOLD_similarity",
			"BOLD_flags",
			"SINTAX_species",
			"SINTAX_bootstrap",
			"Source_of_consensus",
			"Warnings",
		],
	)
	consensus_table.to_csv(consensus_path)
	_annotate_fasta(fasta_path, {row["ASV_ID"]: row for row in consensus_rows}, annotated_fasta)
	files["taxonomy_consensus_csv"] = consensus_path
	files["annotated_fasta"] = annotated_fasta
	files["conflict_report_csv"] = conflict_report
	for source_name, table_path in (("silva", files.get("silva_sintax_csv")), ("pr2", files.get("pr2_sintax_csv"))):
		if table_path:
			log_lines.append(f"{source_name} source output: {table_path}")
	log_lines.append(f"Consensus method: {consensus_method}")
	log_path.write_text("\n".join(log_lines), encoding="utf-8")
	return StepOutput(step_name="taxonomy", success=True, files=files, log=log_path, metrics=metrics)