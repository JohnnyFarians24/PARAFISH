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
from core.preflight import available_memory_bytes, sintax_memory_required

from databases.bold_db import BOLDConfig, run_boldigger as run_boldigger_step
from databases.ncbi_db import BLAST_FAILURE_SIGNATURES, NCBIConfig, run_blast
from databases.pr2_db import PR2Config
from databases.reference_filter import build_marker_reference
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


def _lineage_cache_path() -> Optional[Path]:
	"""A stable, cross-run cache for taxid->lineage lookups, or ``None`` if unwritable.

	Kept beside the reference databases (not inside a run directory) so lineages
	resolved by one run are reused by the next instead of being re-fetched.
	"""
	try:
		from core.config import DB_ROOT

		return Path(DB_ROOT) / "ncbi" / "lineage_cache.json"
	except Exception:
		return None


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


def _is_parquet_path(path: Path) -> bool:
	"""True for ``.parquet`` and BOLDigger3's ``.parquet.snappy`` outputs."""
	name = path.name.lower()
	return name.endswith(".parquet") or name.endswith(".parquet.snappy")


def _read_table_rows(path: Path) -> List[Dict[str, str]]:
	if not path or not path.exists() or path.stat().st_size == 0:
		return []
	if _is_parquet_path(path) or path.suffix.lower() in {".xlsx", ".xls"}:
		try:
			if pd is None:
				raise ImportError("pandas not available")
			if _is_parquet_path(path):
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
	"""Case-insensitive column lookup.

	Source tables disagree on capitalisation: BOLDigger3 emits ``Species``/
	``Genus``/``pct_identity``, BLAST emits ``sscinames``, SINTAX rows are
	normalised to lowercase. An exact/lower/upper probe missed ``Species``
	entirely and silently dropped every BOLD assignment.
	"""
	if key in row:
		return _safe_text(row.get(key))
	target = key.strip().lower().rstrip("_")
	for column, value in row.items():
		if str(column).strip().lower().rstrip("_") == target:
			return _safe_text(value)
	return ""


def _normalize_taxon_value(value: str) -> str:
	clean = value.strip()
	if not clean:
		return ""
	if clean.lower() in {"na", "unclassified", "unknown", "none", "nan", "no-match", "nomatch", "incompletetaxonomy", ""}:
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


# NCBI's scientific names include a lot of placeholder text that is shaped exactly
# like a binomial ("Eukaryota environmental sample") but names no genus at all.
_NON_BINOMIAL_MARKERS = frozenset(
	{"uncultured", "unidentified", "unclassified", "environmental", "sample", "clone", "isolate", "metagenome", "symbiont"}
)


def _genus_from_binomial(species: str) -> str:
	"""``"Micromesistius poutassou"`` -> ``"Micromesistius"``.

	BLAST reports a scientific name and nothing else — no rank breakdown — so an
	otherwise perfectly good species-level hit arrived with an empty genus, and every
	"classified to genus" count and genus-level plot silently ignored it. The genus is
	the first word of a binomial by definition, so recovering it loses no information.
	Anything that is not a real binomial ("uncultured eukaryote", "Bacteria",
	"Eukaryota environmental sample") yields nothing rather than a made-up genus.
	"""
	parts = species.split()
	if len(parts) < 2:
		return ""
	if _NON_BINOMIAL_MARKERS & {part.lower().strip(".") for part in parts}:
		return ""
	genus, epithet = parts[0], parts[1]
	if not genus[:1].isupper() or not genus[1:].islower():
		return ""
	# "sp.", "cf." and "aff." still identify the genus; otherwise the epithet of a
	# real binomial is lowercase.
	if not (epithet.islower() or epithet in {"sp.", "cf.", "aff."}):
		return ""
	return genus


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
		"similarity": _first_nonempty(row, "similarity", "BOLD_similarity", "pct_identity", "pident"),
		# SINTAX is a k-mer classifier: it never aligns, so it has no coverage to
		# report. Only the BLAST-backed sources fill this in.
		"coverage": _first_nonempty(row, "coverage", "qcovhsp", "qcov"),
		"confidence": _first_nonempty(row, "confidence", "SINTAX_bootstrap", "bootstrap", "bitscore"),
		"flags": _first_nonempty(row, "flags", "BOLD_flags", "Warnings"),
		# BLAST's taxid, carried so the NCBI lineage resolver can fill Kingdom..Family.
		"taxid": _first_nonempty(row, "staxids", "taxid"),
	}
	for level in TAXONOMY_LEVELS:
		normalized[level] = _normalize_taxon_value(normalized[level])
	if not normalized["genus"] and normalized["species"]:
		normalized["genus"] = _genus_from_binomial(normalized["species"])
	return normalized


def deepest_rank(record: Dict[str, str]) -> str:
	"""Deepest taxonomic level this record actually assigned, or ``"unassigned"``."""
	for level in reversed(TAXONOMY_LEVELS):
		if _normalize_taxon_value(_safe_text(record.get(level))):
			return level
	return "unassigned"


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


def best_alignment_record(records: list[dict]) -> Optional[dict]:
	"""The strongest alignment-based hit for a feature: highest identity, then coverage.

	SINTAX records are excluded because they carry no identity — they are k-mer
	bootstraps, and ranking them alongside an alignment would compare two different
	things on one scale.

	A record that named no taxon is not a hit, and neither is one at identity 0. BOLD
	writes its non-identifications as rows carrying ``similarity = 0.0`` — which is not
	``None``, so they scored, and with SILVA/PR2 excluded for being SINTAX the *only*
	scoring record was BOLD's non-answer. That is how every feature came to be credited
	``Best_Source=BOLD, Best_Identity=0.0`` while BOLD had identified nothing and the
	species on the row had come from SILVA. The consensus was citing, as its evidence, the
	one source that had no evidence.
	"""
	scored = []
	for record in records:
		identity = _to_float(record.get("similarity"))
		if identity is None or identity <= 0:
			continue
		if not any(record.get(level) for level in TAXONOMY_LEVELS):
			continue
		coverage = _to_float(record.get("coverage"))
		scored.append(((identity, coverage if coverage is not None else -1.0), record))
	if not scored:
		return None
	return max(scored, key=lambda item: item[0])[1]


def best_sintax_record(records: list[dict]) -> Optional[dict]:
	"""The strongest SINTAX (k-mer bootstrap) hit for a feature: highest confidence.

	The mirror of ``best_alignment_record`` for the other evidence family. A SINTAX record
	is a k-mer classifier that NEVER aligns, so it carries a bootstrap confidence and NO
	identity — that absence of identity is what tells it apart from an aligner. The
	distinction matters because the record's ``confidence`` field falls back to an aligner's
	BLAST bitscore (see ``_normalize_source_row``): an NCBI hit therefore carries a
	``confidence`` too, and keying on "has confidence" alone would mislabel its bitscore as a
	SINTAX bootstrap. Requiring identity to be ABSENT excludes NCBI/BOLD cleanly. A record
	that named no taxon is not a hit.
	"""
	scored = []
	for record in records:
		if _to_float(record.get("similarity")) is not None:
			continue  # has an alignment identity → it is an aligner, not SINTAX
		confidence = _to_float(record.get("confidence"))
		if confidence is None:
			continue
		if not any(record.get(level) for level in TAXONOMY_LEVELS):
			continue
		scored.append((confidence, record))
	if not scored:
		return None
	return max(scored, key=lambda item: item[0])[1]


def _build_consensus_row(asv_id: str, records: list[dict], method: str, feature_type: str = "") -> Dict[str, str]:
	method_name = (method or "lca").strip().lower()
	if method_name == "strict":
		consensus = _strict_consensus(records)
	elif method_name == "majority":
		consensus = _majority_consensus(records)
	else:
		consensus = _lca_consensus(records)
	row = {
		"ASV_ID": asv_id,
		"feature_type": feature_type,
		"Consensus_Kingdom": consensus.get("kingdom", ""),
		"Consensus_Phylum": consensus.get("phylum", ""),
		"Consensus_Class": consensus.get("class_", ""),
		"Consensus_Order": consensus.get("order", ""),
		"Consensus_Family": consensus.get("family", ""),
		"Consensus_Genus": consensus.get("genus", ""),
		"Consensus_Species": consensus.get("species", ""),
		"Best_Identity": "",
		"Best_Coverage": "",
		"Best_Source": "",
		# Only sources that actually named a taxon for this feature: a row that reached the
		# table with every rank empty contributed nothing to the LCA and must not be listed
		# as if it had.
		"Consensus_Sources": "|".join(
			sorted({
				str(record.get("source") or "").strip().upper()
				for record in records
				if record.get("source") and any(record.get(level) for level in TAXONOMY_LEVELS)
			})
		),
		"Reliability": "",
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

	# The evidence behind the call. Without these the consensus is an assertion with
	# nothing to judge it by — which is why they now sit in the main table and not
	# only in the per-source file nobody opens.
	consensus_rank = deepest_rank(
		{level: consensus.get(level, "") for level in TAXONOMY_LEVELS}
	)
	best = best_alignment_record(records)
	if best is not None:
		row["Best_Identity"] = _safe_text(best.get("similarity"))
		row["Best_Coverage"] = _safe_text(best.get("coverage"))
		row["Best_Source"] = _safe_text(best.get("source"))
		row["Reliability"] = assignment_reliability(best, consensus_rank)
	else:
		# No aligner hit: grade on the best SINTAX bootstrap instead, and say so.
		sintax = [record for record in records if _to_float(record.get("confidence")) is not None]
		if sintax:
			strongest = max(sintax, key=lambda record: _to_float(record.get("confidence")) or 0.0)
			row["Best_Source"] = _safe_text(strongest.get("source"))
			row["Reliability"] = assignment_reliability(strongest, consensus_rank)
		else:
			row["Reliability"] = "unknown"
	return row


_CONSENSUS_RANK_COLUMNS = [
	"Consensus_Kingdom",
	"Consensus_Phylum",
	"Consensus_Class",
	"Consensus_Order",
	"Consensus_Family",
	"Consensus_Genus",
	"Consensus_Species",
]

# The one definition of the consensus table's shape. Best_Identity / Best_Coverage /
# Reliability sit here rather than only in source_assignments.csv because the
# consensus is the table people actually read, and a taxonomic call without the
# identity and coverage behind it cannot be judged.
CONSENSUS_COLUMNS = [
	"ASV_ID",
	"feature_type",
	*_CONSENSUS_RANK_COLUMNS,
	"Best_Identity",
	"Best_Coverage",
	"Best_Source",
	# Which sources actually fed the LCA for this feature. Distinct from Best_Source,
	# which names only the strongest *alignment* — SINTAX sources carry no identity, so
	# they can never be Best_Source and were therefore invisible in the output despite
	# deciding the call. With this column, "PR2 was harmonised, not dropped" is a fact
	# the table states rather than a claim to be taken on trust.
	"Consensus_Sources",
	"Reliability",
	"BOLD_species",
	"BOLD_similarity",
	"BOLD_flags",
	"SINTAX_species",
	"SINTAX_bootstrap",
	"Source_of_consensus",
	"Warnings",
]


def _annotate_fasta(fasta_path: Path, consensus_rows: Dict[str, Dict[str, str]], output_path: Path) -> None:
	def _write_record(handle, header: str, sequence_lines: list[str]) -> None:
		# Consensus rows are keyed by the bare ID, so the ``;sample=..;size=..``
		# suffixes VSEARCH appends must be stripped before the lookup — keeping
		# them made every annotation fall back to an empty taxonomy.
		asv_id = header.split()[0].split(";")[0]
		taxonomy = consensus_rows.get(asv_id, {})
		size = 1
		if "size=" in header:
			try:
				size = int(header.split("size=", 1)[1].split(";", 1)[0])
			except ValueError:
				size = 1
		tax_parts = [taxonomy.get(key, "") for key in _CONSENSUS_RANK_COLUMNS]
		handle.write(f">{asv_id};size={size};tax={';'.join(tax_parts)}\n")
		for seq_line in sequence_lines:
			handle.write(seq_line + "\n")

	output_path.parent.mkdir(parents=True, exist_ok=True)
	with fasta_path.open("r", encoding="utf-8", errors="replace") as input_handle, output_path.open("w", encoding="utf-8") as output_handle:
		current_header = None
		sequence_lines: list[str] = []
		for raw_line in input_handle:
			line = raw_line.rstrip("\n")
			if line.startswith(">"):
				if current_header is not None:
					_write_record(output_handle, current_header, sequence_lines)
				current_header = line[1:].strip()
				sequence_lines = []
			else:
				sequence_lines.append(line)
		if current_header is not None:
			_write_record(output_handle, current_header, sequence_lines)


# Minimum % identity conventionally required to trust an alignment-based call at
# each rank. Deeper ranks demand more similarity; a 90% hit supports a family, not
# a species. Applied only to sources that actually align (BLAST, BOLD).
RANK_IDENTITY_THRESHOLDS = {
	"species": 97.0,
	"genus": 95.0,
	"family": 90.0,
	"order": 85.0,
	"class_": 80.0,
	"phylum": 75.0,
	"kingdom": 0.0,
}

# SINTAX bootstrap bands. VSEARCH has *already* dropped everything below the run's
# cutoff (0.8 by default), so grading against that same cutoff would mark every
# surviving call "high" — a flag that is true by construction tells the reader
# nothing. These bands sit above it, so the grade varies with the actual bootstrap.
SINTAX_HIGH_CONFIDENCE = 0.95
SINTAX_MODERATE_CONFIDENCE = 0.85

# Minimum SINTAX bootstrap to trust a k-mer call, at any rank (kept for callers).
SINTAX_CONFIDENCE_THRESHOLD = 0.8

# Minimum query coverage for an alignment-based call to count as reliable.
MIN_COVERAGE = 80.0

ASSIGNMENT_COLUMNS = ["ASV_ID", "source", "rank", "taxon", "identity", "coverage", "confidence", "reliability"]


def _to_float(value) -> Optional[float]:
	try:
		return float(str(value).strip())
	except (TypeError, ValueError):
		return None


def assignment_reliability(record: Dict[str, str], rank: str) -> str:
	"""Classify one source's call at ``rank`` as ``high``/``moderate``/``low``.

	Alignment-based sources are judged on identity against a rank-specific
	threshold plus query coverage; SINTAX is judged on its bootstrap. A source
	that reports neither cannot be graded and comes back ``unknown``.
	"""
	if rank == "unassigned":
		return "unassigned"
	identity = _to_float(record.get("similarity"))
	if identity is not None:
		threshold = RANK_IDENTITY_THRESHOLDS.get(rank, 0.0)
		coverage = _to_float(record.get("coverage"))
		if identity < threshold:
			return "low"
		if coverage is not None and coverage < MIN_COVERAGE:
			return "moderate"
		return "high"
	confidence = _to_float(record.get("confidence"))
	if confidence is not None and 0.0 <= confidence <= 1.0:
		if confidence >= SINTAX_HIGH_CONFIDENCE:
			return "high"
		if confidence >= SINTAX_MODERATE_CONFIDENCE:
			return "moderate"
		return "low"
	return "unknown"


def _build_assignment_rows(records_by_asv: Dict[str, List[dict]]) -> List[Dict[str, str]]:
	"""One row per (ASV, source): what it called, how deep, and how much to trust it."""
	rows: List[Dict[str, str]] = []
	for asv_id, records in sorted(records_by_asv.items()):
		for record in records:
			rank = deepest_rank(record)
			rows.append(
				{
					"ASV_ID": asv_id,
					"source": record.get("source", ""),
					"rank": rank,
					"taxon": _safe_text(record.get(rank)) if rank != "unassigned" else "",
					"identity": _safe_text(record.get("similarity")),
					"coverage": _safe_text(record.get("coverage")),
					"confidence": _safe_text(record.get("confidence")),
					"reliability": assignment_reliability(record, rank),
				}
			)
	return rows


def _build_conflict_report(records_by_asv: Dict[str, List[dict]]) -> List[Dict[str, str]]:
	"""Features whose sources disagree at some rank, in a STABLE order.

	Sorted like the consensus is (``sorted(records_by_asv.items())``), and for the same
	reason: this walked the dict in insertion order, which depends on the order the
	sources happened to return records, so two runs over identical data emitted the same
	rows shuffled. A file that changes without its content changing cannot be diffed
	against a reference, which is the only way anyone would notice it changing for a
	reason that matters.
	"""
	report_rows: List[Dict[str, str]] = []
	for asv_id, records in sorted(records_by_asv.items()):
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


def _classify_sintax_source(result_path: Optional[Path], log_path: Optional[Path]) -> tuple[str, str]:
    """Classify a SINTAX source outcome as ``(status, detail)``.

    Distinguishes a genuine "no assignments" from a source that was killed
    (out-of-memory) or otherwise failed — so the UI never renders a silent
    "None" that hides a crashed reference database.
    """
    log_text = ""
    if log_path and Path(log_path).exists():
        try:
            log_text = Path(log_path).read_text(encoding="utf-8", errors="replace")
        except Exception:
            log_text = ""
    lowered = log_text.lower()
    oom = ("terminado pelo sistema" in lowered) or ("sinal 9" in lowered) or ("oom" in lowered) or ("out of memory" in lowered) or ("signal 9" in lowered)
    assigned = 0
    if result_path and Path(result_path).exists() and Path(result_path).stat().st_size > 0:
        try:
            with Path(result_path).open("r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.strip() and not line.lower().startswith("asv_id"):
                        assigned += 1
        except Exception:
            assigned = 0
    if assigned > 0:
        return "ok", f"{assigned} sequências classificadas"
    if oom:
        return "failed", "sem memória (OOM) — base demasiado grande para a RAM disponível"
    if "reference preparation failed" in lowered or "error" in lowered:
        return "failed", "erro ao preparar/correr a base — ver log"
    return "no_hits", "correu mas sem atribuições taxonómicas"


def _count_fasta_records(fasta_path) -> int:
	"""How many features were submitted to the sources — the denominator in a status."""
	try:
		with Path(fasta_path).open("r", errors="replace") as handle:
			return sum(1 for line in handle if line.startswith(">"))
	except OSError:
		return 0


def _source_assignment_counts(records_by_asv) -> Dict[str, int]:
	"""Features each source actually assigned a taxon to, keyed by record source name.

	``records_by_asv`` is the single source of truth here: it is what
	``source_assignments.csv`` and the consensus are built from. A record only counts
	if it carries a real taxon at some rank — a row that reached the table with every
	rank empty (BOLD's "unassigned, identity 0.0") is not an assignment.
	"""
	assigned: Dict[str, set] = defaultdict(set)
	for asv_id, records in (records_by_asv or {}).items():
		for record in records:
			source = str(record.get("source") or "").strip().upper()
			if source and any(record.get(level) for level in TAXONOMY_LEVELS):
				assigned[source].add(asv_id)
	return {source: len(ids) for source, ids in assigned.items()}


# The markers each source's reference can actually classify. Only sources with a real
# restriction appear here; absence means "no known restriction", not "anything".
#
# BOLD: ``bold_params["db"] = 1`` is BOLD's public *animal* library, which is COI. Asked
# to identify 18S amplicons it returns rows with every rank empty — and "returned nothing"
# is the same shape as an honest absence of matches, which is how it came to be reported
# as an ordinary ``no_hits`` on data it was never able to answer for.
# SILVA/PR2: both shipped releases are SSU (``SILVA_138.2_SSURef_NR99``,
# ``pr2_5.1.1_SSU``) — small-subunit rRNA. They cannot carry a protein-coding COI
# amplicon any more than BOLD can carry 18S, but only BOLD's direction was ever
# encoded, so a COI run put SILVA and PR2 through the amplicon cut and the SINTAX
# load and then reported the nothing they found as an ordinary ``no_hits``. PR2 is
# the eukaryotic SSU, so 16S runs are outside it too.
SOURCE_SUPPORTED_MARKERS = {
	"BOLD": {"COI"},
	"SILVA": {"16S", "18S"},
	"PR2": {"18S"},
}


def _log_failure_signature(log_path) -> str:
	"""The BLAST failure phrase present in ``log_path``, or ``""``.

	Reuses ``ncbi_db.BLAST_FAILURE_SIGNATURES`` rather than restating them, because the
	tool being judged is the same blastn: it exits 0 after printing "(Severe Error) ...
	search failed", so the exit code says nothing and the log says everything.
	"""
	if not log_path:
		return ""
	try:
		text = Path(log_path).read_text(encoding="utf-8", errors="replace")
	except OSError:
		return ""
	lowered = text.lower()
	for signature in BLAST_FAILURE_SIGNATURES:
		if str(signature).lower() in lowered:
			return str(signature)
	return ""


def _marker_incompatible(source_key: str, marker) -> bool:
	"""True when ``source_key``'s reference cannot carry ``marker``.

	Unknown markers are never incompatible: "we could not tell what this is" is not
	evidence that a source is inapplicable, and refusing to run on it would turn a
	gap in configuration into a silent loss of a database.
	"""
	supported = SOURCE_SUPPORTED_MARKERS.get(str(source_key or "").upper())
	marker_upper = str(marker or "").strip().upper()
	if not supported or not marker_upper or marker_upper == "OTHER":
		return False
	return marker_upper not in supported


def _not_applicable(source: str, marker) -> Dict[str, str]:
	"""The status entry for a source skipped before it ever ran.

	No ``key``/``raw``: those are what ``_reconcile_source_status`` judges a source by,
	and there is nothing to judge here. Without a key the entry passes through
	reconciliation untouched, which is what keeps this verdict final.
	"""
	supported = SOURCE_SUPPORTED_MARKERS.get(str(source or "").upper()) or set()
	return {
		"source": source,
		"status": "not_applicable",
		"detail": (
			f"{source} indexa {'/'.join(sorted(supported))}; este run é "
			f"{str(marker).strip().upper()} — a fonte foi saltada sem ser filtrada nem carregada"
		),
	}

# Seed for VSEARCH's PRNG in SINTAX. Any fixed non-zero value works; what matters is that
# it is not 0, which is vsearch's default and means "seed from a random source". A
# classification that changes between two runs of identical data is not a result, it is a
# sample from a distribution — and a reviewer who re-runs the pipeline must get our numbers.
SINTAX_RANDSEED = 42


def _reconcile_source_status(
	source_status: List[Dict[str, str]], records_by_asv, feature_total: int, marker: str = ""
) -> List[Dict[str, str]]:
	"""Force every source's status to agree with what it actually produced.

	Each source used to judge itself by a different proxy: BOLD by whether its result
	*file* existed, NCBI by BLAST's ``hit_count``, custom BLAST by a match counter.
	None of those answers the question the column claims to answer — "did this source
	assign a taxon to anything" — so BOLD reported ``ok, identificação concluída`` on a
	run where all 7 of its rows were unassigned with identity 0.0, and custom BLAST
	reported ``no_hits`` off a raw output of 0 bytes, which is a failure, not an
	absence of matches.

	The rule: no source claims success without at least one assignment. With none:

	* the source could not carry this run's marker  -> ``not_applicable``
	* empty/absent raw output                       -> ``failed``
	* real output                                   -> ``no_hits``

	``not_applicable`` exists because ``no_hits`` lies by omission on a marker mismatch.
	"BOLD found nothing in these 18S reads" and "BOLD does not index 18S at all" are the
	same word but different facts: the first invites the reader to conclude the taxa are
	absent from BOLD, when BOLD was never asked a question it could answer. Distinguishing
	them is the whole point of this column.
	"""
	counts = _source_assignment_counts(records_by_asv)
	marker_upper = str(marker or "").strip().upper()
	reconciled: List[Dict[str, str]] = []
	for entry in source_status:
		key = str(entry.pop("key", "") or "").upper()
		raw = entry.pop("raw", None)
		log = entry.pop("log", None)
		status = entry.get("status")
		# Judge every source that ran, not only the ones claiming "ok". NCBI reports
		# "ok if hits > 0 else no_hits" up front, so a remote BLAST that died outright
		# ("CPU usage limit was exceeded" — 0-byte result, error-filled log) arrived here
		# already labelled no_hits and was waved through: a failure wearing the words of
		# an honest absence. Only "skipped" is left alone, because a source that never ran
		# has nothing to reconcile against.
		if key and status != "skipped":
			assigned = counts.get(key, 0)
			supported = SOURCE_SUPPORTED_MARKERS.get(key)
			incompatible = bool(supported and marker_upper and marker_upper not in supported)
			# With BLAST, zero hits legitimately produces a ZERO-BYTE table: the search ran,
			# matched nothing, and wrote nothing. So an empty file cannot mean "failed" on
			# its own — custom BLAST reported `failed` on a run whose besthit CSV was
			# complete and correct, every feature honestly marked "None". The evidence for
			# failure is in the LOG, where blastn prints "(Severe Error) ... search failed"
			# while still exiting 0; that is the same signal ncbi_db already judges by.
			# File size is used only as a fallback, for sources that ship no log.
			failure = _log_failure_signature(log)
			produced_output = (
				not failure
				if log is not None
				else bool(raw and Path(raw).exists() and Path(raw).stat().st_size > 0)
			)
			if assigned:
				entry["status"] = "ok"
				entry["detail"] = f"{assigned}/{feature_total} features com taxonomia · {entry.get('detail', '')}".strip(" ·")
			elif incompatible:
				entry["status"] = "not_applicable"
				entry["detail"] = (
					f"{key} indexa {'/'.join(sorted(supported))}; este run é {marker_upper} — "
					f"a fonte não podia identificar estas features, e 0/{feature_total} não é uma ausência real"
				)
			elif produced_output:
				entry["status"] = "no_hits"
				entry["detail"] = f"correu mas não atribuiu taxonomia a nenhuma das {feature_total} features"
			else:
				entry["status"] = "failed"
				entry["detail"] = (
					f"a ferramenta reportou um erro: {failure}"
					if failure
					else "não produziu resultado utilizável (saída vazia ou ausente)"
				)
		reconciled.append(entry)
	return reconciled


def _write_csv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
	path.parent.mkdir(parents=True, exist_ok=True)
	with path.open("w", newline="", encoding="utf-8") as handle:
		writer = csv.DictWriter(handle, fieldnames=columns)
		writer.writeheader()
		writer.writerows(rows)


def _sintax_memory_shortfall(db_path: Path) -> Optional[str]:
	"""Message explaining why this reference cannot be loaded, or ``None`` if it can.

	VSEARCH --sintax holds the reference sequences AND a k-mer index in RAM, so the
	requirement is a multiple of the on-disk FASTA, not its size. Comparing against
	the raw file size (as this used to) let a 786 MB SILVA through on a box with
	1.0 GB free — and it was then SIGKILLed partway through, which reads to the user
	as an unexplained failure. Refusing up front is the same outcome, minus the wait
	and the confusion.
	"""
	needed = sintax_memory_required(Path(db_path))
	available = available_memory_bytes()
	if needed is None or available is None or needed <= available:
		return None
	size = Path(db_path).stat().st_size
	return (
		f"ignorada: a base ocupa {size / 1e9:.1f} GB e o SINTAX precisa de cerca de "
		f"{needed / 1e9:.1f} GB de RAM, mas so ha {available / 1e9:.1f} GB disponiveis. "
		f"Usa uma referencia mais pequena (filtrada para o teu marcador), ou mais RAM."
	)


def _apply_marker_filter(
	db_path: Path,
	marker,
	primers: tuple,
	source_name: str,
	log_lines: List[str],
	threads: int = 2,
) -> tuple[Path, Optional[str]]:
	"""Swap in a marker-sized copy of ``db_path``. Returns ``(path, note)``.

	This runs *before* the memory gate on purpose: the gate's job is to answer "does
	this reference fit", and the reference it should be asked about is the one that is
	actually going to be loaded. Filtering first is what turns SILVA from a 2.5 GB
	skip into a 0.1 GB run.

	Failure here is never fatal — the original path comes back and the gate decides, as
	it did before this existed.
	"""
	try:
		reference = build_marker_reference(Path(db_path), marker, primers, threads=threads)
	except Exception as exc:  # a reference that cannot be shrunk is not a failed run
		log_lines.append(f"{source_name}: falha ao filtrar a referência ({exc}); a usar a base completa.")
		return Path(db_path), None
	if reference is None:
		return Path(db_path), None
	note = f"{source_name}: {reference.summary()}"
	log_lines.append(note)
	return reference.path, reference.summary()


def run_sintax(fasta_path: Path, db_path: Path, cutoff: float, threads: int, output_dir: Optional[Path] = None, prefix: str = "sintax") -> Path:
	"""Run VSEARCH SINTAX and return the result CSV path."""
	output_dir = output_dir or (fasta_path.parent / "taxonomy")
	output_dir.mkdir(parents=True, exist_ok=True)
	result_path = output_dir / f"{prefix}_result.csv"
	log_path = output_dir / f"{prefix}.log"

	preamble = _sintax_memory_shortfall(Path(db_path))

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
		# SINTAX bootstraps by randomly subsampling k-mers, and vsearch's --randseed
		# defaults to 0 = "use a random data source". Unseeded, the same reads against the
		# same reference give different bootstraps on every run: SILVA called OTU4 at 0.99
		# and then 0.98 on a re-run of identical input. Those numbers decide `Reliability`
		# and can cross `--sintax_cutoff`, so an unseeded run means a taxon can appear or
		# vanish between two runs of the same data — and nobody could reproduce a result.
		"--randseed",
		str(SINTAX_RANDSEED),
		"--log",
		str(log_path),
	]
	try:
		completed = subprocess.run(command, capture_output=True, text=True, check=True)
		body = (completed.stdout or "") + ("\n" if completed.stdout and completed.stderr else "") + (completed.stderr or "")
		log_path.write_text(((preamble + "\n\n") if preamble else "") + body, encoding="utf-8")
	except subprocess.CalledProcessError as exc:
		if exc.returncode is not None and exc.returncode < 0:
			reason = (
				f"VSEARCH SINTAX foi terminado pelo sistema (sinal {-exc.returncode}), "
				f"quase de certeza por falta de memória (OOM)."
			)
			message = (preamble + "\n\n" if preamble else "") + reason
		else:
			message = (preamble + "\n\n" if preamble else "") + str(exc)
		log_path.write_text(message, encoding="utf-8")
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
	return ConsensusTable(consensus_rows, list(CONSENSUS_COLUMNS))


_SINTAX_RANK_BY_PREFIX = {
	"d": "kingdom",
	"k": "kingdom",
	"p": "phylum",
	"c": "class_",
	"o": "order",
	"f": "family",
	"g": "genus",
	"s": "species",
}


def _split_taxon_confidence(value: str) -> tuple[str, str]:
	"""Split ``Genus_species(0.98)`` into ``("Genus_species", "0.98")``.

	Only a *numeric* parenthetical is a bootstrap. This used to strip whatever came
	last in parentheses, which is wrong twice over on SILVA: 7.7% of its eukaryote
	names carry a common name (``Gadus_morhua_(Atlantic_cod)``), and the cutoff-filtered
	column carries no bootstraps at all — so that column's species were silently
	renamed to ``Gadus_morhua_`` and their bootstrap became the word ``Atlantic_cod``.
	"""
	value = value.strip()
	if not (value.endswith(")") and "(" in value):
		return value, ""
	head, _, tail = value.rpartition("(")
	candidate = tail.rstrip(")").strip()
	try:
		float(candidate)
	except ValueError:
		return value, ""
	return head.strip(), candidate


def _parse_sintax_assignment(field: str) -> tuple[Dict[str, str], str]:
	"""Parse a VSEARCH SINTAX taxonomy field into ranks plus the best bootstrap.

	The field looks like ``d:Eukaryota(1.00),p:Obazoa(1.00),...,s:Teleostei(0.64)``.
	Bootstraps are optional (the cutoff-filtered column has none).
	"""
	ranks: Dict[str, str] = {}
	last_confidence = ""
	for chunk in field.split(","):
		chunk = chunk.strip()
		if not chunk or ":" not in chunk:
			continue
		prefix, _, value = chunk.partition(":")
		rank = _SINTAX_RANK_BY_PREFIX.get(prefix.strip().lower())
		if not rank:
			continue
		value, confidence = _split_taxon_confidence(value)
		value = _normalize_taxon_value(value)
		if not value:
			continue
		ranks[rank] = value
		if confidence:
			last_confidence = confidence
	return ranks, last_confidence


def _sintax_confidence_for(raw_field: str, rank: str) -> str:
	"""Return the bootstrap SINTAX attached to ``rank`` in its raw prediction field."""
	for chunk in raw_field.split(","):
		chunk = chunk.strip()
		prefix, _, value = chunk.partition(":")
		if _SINTAX_RANK_BY_PREFIX.get(prefix.strip().lower()) != rank:
			continue
		# A name ending in "(Atlantic_cod)" is not a bootstrap of "Atlantic_cod".
		return _split_taxon_confidence(value)[1]
	return ""


def _collect_sintax_records(source: str, table_path: Optional[Path]) -> Dict[str, List[dict]]:
	"""Read a VSEARCH ``--tabbedout`` SINTAX table into normalized records.

	VSEARCH writes a headerless TSV (query, all predictions with bootstraps,
	strand, predictions passing the cutoff). Feeding it to ``csv.DictReader``
	consumed the first hit as a header row and yielded no usable ranks, so no
	SINTAX source ever reached the consensus.
	"""
	rows_by_asv: Dict[str, List[dict]] = defaultdict(list)
	if not (table_path and table_path.exists() and table_path.stat().st_size > 0):
		return rows_by_asv
	with table_path.open("r", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			line = line.rstrip("\n")
			if not line.strip() or line.lower().startswith("asv_id"):
				continue
			fields = line.split("\t")
			asv_id = fields[0].split(";")[0].strip()
			if not asv_id:
				continue
			# Prefer the cutoff-filtered column (4th) and fall back to the raw
			# predictions (2nd) so a row is never dropped for lack of bootstraps.
			filtered = fields[3] if len(fields) > 3 else ""
			raw = fields[1] if len(fields) > 1 else ""
			ranks, confidence = _parse_sintax_assignment(filtered or raw)
			if not ranks:
				continue
			if not confidence and raw:
				# Only the raw column carries bootstraps; report the one attached to
				# the deepest rank that survived the cutoff.
				deepest = max(ranks, key=TAXONOMY_LEVELS.index)
				confidence = _sintax_confidence_for(raw, deepest)
			# SINTAX reports a bootstrap only: no alignment means no identity and no
			# coverage, and inventing either would misrepresent the evidence.
			record = {"asv_id": asv_id, "source": source, "confidence": confidence, "similarity": "", "coverage": "", "flags": ""}
			for level in TAXONOMY_LEVELS:
				record[level] = ranks.get(level, "")
			rows_by_asv[asv_id].append(record)
	return rows_by_asv


def _collect_blast_hits(source: str, table_path: Optional[Path], outfmt: str) -> Dict[str, List[dict]]:
	"""Read a headerless ``blastn -outfmt 6 ...`` table into ALL normalized hits.

	The column names come from the outfmt string used to produce the file, since
	BLAST never writes a header line.
	"""
	hits_by_asv: Dict[str, List[dict]] = defaultdict(list)
	if not (table_path and table_path.exists() and table_path.stat().st_size > 0):
		return hits_by_asv
	columns = [token for token in str(outfmt or "").split() if not token.isdigit()]
	if not columns:
		return hits_by_asv
	with table_path.open("r", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			line = line.rstrip("\n")
			if not line.strip() or line.startswith("#"):
				continue
			fields = line.split("\t")
			row = dict(zip(columns, fields))
			normalized = _normalize_source_row(source, row)
			if not normalized["asv_id"]:
				continue
			normalized["subject"] = _first_nonempty(row, "sseqid", "saccver")
			normalized["evalue"] = _first_nonempty(row, "evalue")
			normalized["bitscore"] = _first_nonempty(row, "bitscore")
			hits_by_asv[normalized["asv_id"]].append(normalized)
	return hits_by_asv


def _collect_blast_records(source: str, table_path: Optional[Path], outfmt: str) -> Dict[str, List[dict]]:
	"""Best BLAST hit per query, for the consensus.

	BLAST reports several hits per query; only the strongest reaches the consensus so
	that one query does not outvote the other sources. The discarded hits are not
	lost — ``_collect_blast_hits`` keeps them for the ranked hits table, which is what
	lets a user pick the second-best hit when it is the more plausible one.
	"""
	rows_by_asv: Dict[str, List[dict]] = defaultdict(list)
	for asv_id, hits in _collect_blast_hits(source, table_path, outfmt).items():
		best = max(hits, key=lambda record: _to_float(record.get("bitscore")) or 0.0)
		rows_by_asv[asv_id].append(best)
	return rows_by_asv


HITS_COLUMNS = [
	"ASV_ID",
	"feature_type",
	"hit_rank",
	"source",
	"rank",
	"taxon",
	"species",
	"genus",
	"family",
	"identity",
	"coverage",
	"confidence",
	"reliability",
	"subject",
	"evalue",
	"bitscore",
]


def _build_hits_rows(hits_by_asv: Dict[str, List[dict]], feature_type: str = "") -> List[Dict[str, str]]:
	"""Every candidate hit for every feature, ranked best-first within each feature.

	Ranking is by % identity then query coverage, so ``hit_rank == 1`` is the best
	supported call. The runners-up are kept deliberately: the top BLAST hit is often
	an unannotated environmental sequence while hit 2 or 3 is the named species the
	user actually needs, and that choice belongs to them, not to us.
	"""
	rows: List[Dict[str, str]] = []
	for asv_id, hits in sorted(hits_by_asv.items()):
		def _sort_key(record: dict):
			identity = _to_float(record.get("similarity"))
			coverage = _to_float(record.get("coverage"))
			confidence = _to_float(record.get("confidence"))
			return (
				identity if identity is not None else -1.0,
				coverage if coverage is not None else -1.0,
				confidence if confidence is not None else -1.0,
			)

		for index, record in enumerate(sorted(hits, key=_sort_key, reverse=True), start=1):
			rank = deepest_rank(record)
			rows.append(
				{
					"ASV_ID": asv_id,
					"feature_type": feature_type,
					"hit_rank": index,
					"source": record.get("source", ""),
					"rank": rank,
					"taxon": _safe_text(record.get(rank)) if rank != "unassigned" else "",
					"species": _safe_text(record.get("species")),
					"genus": _safe_text(record.get("genus")),
					"family": _safe_text(record.get("family")),
					"identity": _safe_text(record.get("similarity")),
					"coverage": _safe_text(record.get("coverage")),
					"confidence": _safe_text(record.get("confidence")),
					"reliability": assignment_reliability(record, rank),
					"subject": _safe_text(record.get("subject")),
					"evalue": _safe_text(record.get("evalue")),
					"bitscore": _safe_text(record.get("bitscore")),
				}
			)
	return rows


# One row per feature: the call, and the two evidence families that back it, side by side.
# NO cross-metric winner — alignment identity and SINTAX bootstrap live in SEPARATE columns
# so a reader sees whether the two families agree, and is never handed a single blended score
# that hides a disagreement (the Phase-3 bug). This is the supervisors' deliverable table.
SUMMARY_COLUMNS = [
	"ASV_ID",
	"feature_type",
	"consensus_taxon",
	"consensus_rank",
	"consensus_reliability",
	"align_source",
	"align_taxon",
	"align_identity",
	"sintax_source",
	"sintax_taxon",
	"sintax_confidence",
]

# Consensus columns from deepest to shallowest — the summary's single taxon call is the
# deepest rank the LCA actually resolved, read from the consensus row already built (one
# source), never recomputed here.
_SUMMARY_CONSENSUS_RANKS = [
	("Consensus_Species", "species"),
	("Consensus_Genus", "genus"),
	("Consensus_Family", "family"),
	("Consensus_Order", "order"),
	("Consensus_Class", "class"),
	("Consensus_Phylum", "phylum"),
	("Consensus_Kingdom", "kingdom"),
]


def _consensus_call(consensus_row: Dict[str, str]) -> tuple:
	"""The deepest resolved rank of the consensus row, as ``(taxon, rank)``."""
	for column, rank in _SUMMARY_CONSENSUS_RANKS:
		value = _safe_text(consensus_row.get(column))
		if value:
			return value, rank
	return "", "unassigned"


def _record_taxon(record: Optional[dict]) -> str:
	"""The named taxon of a hit at its deepest resolved rank — same rule the hits table uses."""
	if record is None:
		return ""
	rank = deepest_rank(record)
	return _safe_text(record.get(rank)) if rank != "unassigned" else ""


def _build_summary_rows(records_by_asv: Dict[str, List[dict]], consensus_rows: List[Dict[str, str]], feature_type: str = "") -> List[Dict[str, str]]:
	"""One row per feature: consensus call plus the best hit of EACH evidence family.

	Derived from the same in-memory records and consensus rows the other taxonomy tables
	are — this file is written by the pipeline and rendered as-is by the UI, so the screen
	and the download cannot diverge.
	"""
	by_id = {row.get("ASV_ID"): row for row in consensus_rows}
	rows: List[Dict[str, str]] = []
	for asv_id, records in sorted(records_by_asv.items()):
		consensus_row = by_id.get(asv_id, {})
		taxon, rank = _consensus_call(consensus_row)
		align = best_alignment_record(records)
		sintax = best_sintax_record(records)
		rows.append(
			{
				"ASV_ID": asv_id,
				"feature_type": feature_type,
				"consensus_taxon": taxon,
				"consensus_rank": rank,
				"consensus_reliability": _safe_text(consensus_row.get("Reliability")),
				"align_source": _safe_text(align.get("source")) if align else "",
				"align_taxon": _record_taxon(align),
				"align_identity": _safe_text(align.get("similarity")) if align else "",
				"sintax_source": _safe_text(sintax.get("source")) if sintax else "",
				"sintax_taxon": _record_taxon(sintax),
				"sintax_confidence": _safe_text(sintax.get("confidence")) if sintax else "",
			}
		)
	return rows


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
	custom_db_kind: Optional[str] = None,
	abundance_matrix_path: Optional[Path] = None,
	consensus_method: str = "lca",
	threads: int = 8,
	sintax_cutoff: float = 0.8,
	feature_set: str = "otu",
	marker: str = "",
	primers: tuple = ("", ""),
) -> StepOutput:
	"""Run all active taxonomy sources for ONE feature set and build its consensus.

	``marker`` and ``primers`` come from the run's own configuration and are used to cut
	SILVA/PR2 down to the region this run amplified (see ``databases.reference_filter``).
	Left empty, the references are used whole — which on a small machine means the memory
	gate skips them.

	``feature_set`` is ``"otu"`` or ``"asv"``. The caller invokes this once per
	feature set, with a distinct ``output_dir`` each time, so OTU and ASV taxonomies
	are produced by the identical code path and land in separate files. They are
	never combined into one table.
	"""
	output_dir.mkdir(parents=True, exist_ok=True)
	feature_type = "ASV" if str(feature_set).lower() == "asv" else "OTU"
	log_path = output_dir / "taxonomy_pipeline.log"
	consensus_path = output_dir / "taxonomy_consensus.csv"
	annotated_fasta = output_dir / "annotated.fasta"
	conflict_report = output_dir / "conflict_report.csv"
	hits_path = output_dir / "taxonomy_hits.csv"
	files: Dict[str, Path] = {"log": log_path}
	metrics: Dict[str, float] = {}
	log_lines: List[str] = []
	records_by_asv: Dict[str, List[dict]] = defaultdict(list)
	# Every candidate hit, not just the winner: this is what backs the ranked hits
	# table and the "give me the second-best" workflow.
	hits_by_asv: Dict[str, List[dict]] = defaultdict(list)
	active_sources = 0
	source_status: List[Dict[str, str]] = []

	def _merge_rows(rows: Dict[str, List[dict]], hits: Optional[Dict[str, List[dict]]] = None) -> None:
		nonlocal active_sources
		if rows:
			active_sources += 1
		for asv_id, source_rows in rows.items():
			records_by_asv[asv_id].extend(source_rows)
		# A source with no separate hit list contributes its consensus records as
		# its hits (SINTAX yields exactly one call per feature).
		for asv_id, hit_rows in (hits if hits is not None else rows).items():
			hits_by_asv[asv_id].extend(hit_rows)

	def _merge_records(source_name: str, table_path: Optional[Path]) -> None:
		_merge_rows(_collect_source_records(source_name, table_path))

	def _merge_sintax(source_name: str, table_path: Optional[Path]) -> None:
		_merge_rows(_collect_sintax_records(source_name, table_path))

	if silva_config is not None and _marker_incompatible("SILVA", marker):
		# Before any of the work: sintax_ready_path() alone rewrites and caches a
		# multi-hundred-MB reference, and the amplicon cut after it costs minutes.
		entry = _not_applicable("SILVA", marker)
		metrics["silva_not_applicable"] = 1.0
		log_lines.append(entry["detail"])
		source_status.append(entry)
	elif silva_config is not None:
		try:
			silva_db_path = silva_config.sintax_ready_path()
		except Exception as exc:
			metrics["silva_error"] = 1.0
			log_lines.append(f"SILVA reference preparation failed: {exc}")
			source_status.append({"source": "SILVA", "status": "failed", "detail": f"erro a preparar a base: {exc}"})
		else:
			# Cut the reference to this run's marker before asking whether it fits: a
			# full SILVA never fits on a small box, but its 18S V4 slice comfortably does.
			silva_db_path, filter_note = _apply_marker_filter(silva_db_path, marker, primers, "SILVA", log_lines, threads=threads)
			# VSEARCH --sintax loads the whole reference into RAM. Running it anyway on a
			# machine that cannot hold it means the OS kills it mid-run (SIGKILL) after
			# burning minutes; skipping up front is faster and says why.
			shortfall = _sintax_memory_shortfall(silva_db_path)
			if shortfall:
				metrics["silva_skipped"] = 1.0
				log_lines.append(shortfall)
				source_status.append({"source": "SILVA", "status": "skipped", "detail": shortfall})
			else:
				silva_path = run_sintax(fasta_path, silva_db_path, sintax_cutoff, threads, output_dir=output_dir / "silva", prefix="silva_sintax")
				files["silva_sintax_csv"] = silva_path
				_merge_sintax("SILVA", silva_path)
				status, detail = _classify_sintax_source(silva_path, output_dir / "silva" / "silva_sintax.log")
				if filter_note:
					detail = f"{detail} · {filter_note}"
				source_status.append({"source": "SILVA", "status": status, "detail": detail, "key": "SILVA", "raw": silva_path})

	if pr2_config is not None and _marker_incompatible("PR2", marker):
		entry = _not_applicable("PR2", marker)
		metrics["pr2_not_applicable"] = 1.0
		log_lines.append(entry["detail"])
		source_status.append(entry)
	elif pr2_config is not None:
		try:
			pr2_db_path = pr2_config.sintax_ready_path()
		except Exception as exc:
			metrics["pr2_error"] = 1.0
			log_lines.append(f"PR2 reference preparation failed: {exc}")
			source_status.append({"source": "PR2", "status": "failed", "detail": f"erro a preparar a base: {exc}"})
		else:
			# PR2 is eukaryote-only already, so here it is the amplicon cut, not the
			# lineage filter, that does the shrinking.
			pr2_db_path, filter_note = _apply_marker_filter(pr2_db_path, marker, primers, "PR2", log_lines, threads=threads)
			# Same memory gate as SILVA: PR2 had none, so on a small machine it was
			# started and then SIGKILLed, which looked like a mysterious failure.
			shortfall = _sintax_memory_shortfall(pr2_db_path)
			if shortfall:
				metrics["pr2_skipped"] = 1.0
				log_lines.append(shortfall)
				source_status.append({"source": "PR2", "status": "skipped", "detail": shortfall})
			else:
				pr2_path = run_sintax(fasta_path, pr2_db_path, sintax_cutoff, threads, output_dir=output_dir / "pr2", prefix="pr2_sintax")
				files["pr2_sintax_csv"] = pr2_path
				_merge_sintax("PR2", pr2_path)
				status, detail = _classify_sintax_source(pr2_path, output_dir / "pr2" / "pr2_sintax.log")
				if filter_note:
					detail = f"{detail} · {filter_note}"
				source_status.append({"source": "PR2", "status": status, "detail": detail, "key": "PR2", "raw": pr2_path})

	custom_blast_success = False
	if custom_db_path is not None and custom_db_path.exists():
		kind = custom_db_kind
		if not kind:
			from databases.custom_db import detect_fasta_kind

			kind, _ = detect_fasta_kind(custom_db_path)
		if kind == "blast":
			# Plain reference FASTA: build a BLAST db and align the OTU/ASV
			# representatives against it, then cross with the abundance matrix.
			from analysis.custom_db_compare import build_custom_comparison, write_detection_table
			from databases.custom_db import run_custom_blast

			custom_result = run_custom_blast(
				fasta_path,
				custom_db_path,
				output_dir / "custom_blast",
				perc_identity=float(90.0),
				threads=threads,
			)
			if custom_result.success:
				custom_blast_success = True
				files.update(custom_result.files)
				metrics["custom_matched"] = custom_result.metrics.get("matched", 0.0)
				metrics["custom_total"] = custom_result.metrics.get("total", 0.0)
				besthit = custom_result.files.get("custom_blast_besthit_csv")
				if abundance_matrix_path and besthit and Path(abundance_matrix_path).exists():
					comparison = build_custom_comparison(besthit, abundance_matrix_path)
					if not comparison.get("error"):
						detection_csv = write_detection_table(comparison.get("detection") or {}, output_dir / "custom_blast" / "custom_detection.csv")
						if detection_csv:
							files["custom_detection_csv"] = detection_csv
				log_lines.append(f"Custom BLAST comparison: {int(metrics.get('custom_matched', 0))}/{int(metrics.get('custom_total', 0))} OTUs com match.")
				matched_n = int(metrics.get("custom_matched", 0))
				total_n = int(metrics.get("custom_total", 0))
				# A 0-byte custom_blast_raw.tsv reported "no_hits" — but no output at
				# all is a broken database or a BLAST that never ran, not an honest
				# absence of matches. Reconciliation tells the two apart.
				source_status.append({
					"source": "Base customizada (BLAST)",
					"status": "ok" if matched_n > 0 else "no_hits",
					"detail": f"{matched_n}/{total_n} OTUs com correspondência",
					"key": "CUSTOM",
					"raw": custom_result.files.get("custom_blast_raw_tsv"),
					# The log, not the raw table's size, decides failure vs absence: a
					# clean blastn that matched nothing writes a 0-byte table, which the
					# size heuristic called "failed" on a run whose besthit CSV was
					# complete with every feature honestly marked "None".
					"log": custom_result.files.get("log"),
				})
			else:
				metrics["custom_error"] = 1.0
				log_lines.append(custom_result.error_message or "Custom BLAST falhou")
				source_status.append({"source": "Base customizada (BLAST)", "status": "failed", "detail": custom_result.error_message or "falhou"})
		else:
			custom_path = run_sintax(fasta_path, custom_db_path, sintax_cutoff, threads, output_dir=output_dir / "custom", prefix="custom_sintax")
			files["custom_sintax_csv"] = custom_path
			_merge_sintax("CUSTOM", custom_path)
			status, detail = _classify_sintax_source(custom_path, output_dir / "custom" / "custom_sintax.log")
			source_status.append({"source": "Base customizada (SINTAX)", "status": status, "detail": detail, "key": "CUSTOM", "raw": custom_path})

	if ncbi_config is not None:
		ncbi_result = run_blast(fasta_path, ncbi_config, output_dir / "ncbi")
		if ncbi_result.success:
			files.update(ncbi_result.files)
			blast_tsv = ncbi_result.files.get("blast_tsv")
			if isinstance(blast_tsv, Path) and blast_tsv.exists():
				ncbi_records = _collect_blast_records("NCBI", blast_tsv, ncbi_config.outfmt)
				ncbi_hits = _collect_blast_hits("NCBI", blast_tsv, ncbi_config.outfmt)
				# BLAST gives only the binomial; resolve the rest of the lineage
				# (Kingdom..Family) from each hit's taxid. Best-effort and cached.
				# Both the consensus records and the ranked hits go into the flat
				# list so a single lineage lookup enriches every dict in both.
				flat_ncbi = [rec for recs in ncbi_records.values() for rec in recs]
				flat_ncbi += [rec for recs in ncbi_hits.values() for rec in recs]
				taxids = [rec.get("taxid") for rec in flat_ncbi if rec.get("taxid")]
				if taxids:
					from databases.ncbi_lineage import fill_records, resolve_lineages

					lineages, lineage_warnings = resolve_lineages(
						taxids, cache_path=_lineage_cache_path()
					)
					fill_records(flat_ncbi, lineages)
					log_lines.extend(lineage_warnings)
				_merge_rows(ncbi_records, ncbi_hits)
				hits = int((getattr(ncbi_result, "metrics", {}) or {}).get("hit_count", 0))
				# hit_count counts BLAST rows; a hit with no resolvable species name is not
				# an assignment. Reconciliation decides.
				source_status.append({"source": "NCBI", "status": "ok" if hits > 0 else "no_hits", "detail": f"{hits} hits BLAST", "key": "NCBI", "raw": blast_tsv})
			else:
				metrics["ncbi_error"] = 1.0
				log_lines.append("NCBI result missing blast_tsv path.")
				source_status.append({"source": "NCBI", "status": "failed", "detail": "resultado BLAST em falta"})
		else:
			metrics["ncbi_error"] = 1.0
			log_lines.append(ncbi_result.error_message or "NCBI failed")
			source_status.append({"source": "NCBI", "status": "failed", "detail": ncbi_result.error_message or "falhou"})

	if bold_config is not None and _marker_incompatible("BOLD", marker):
		# BOLD's incompatibility was already known, but only after the fact: the run
		# still paid for a full remote identification before the reconciliation
		# relabelled the empty result. Same rule as SILVA/PR2 — decide before the work.
		entry = _not_applicable("BOLD", marker)
		metrics["bold_not_applicable"] = 1.0
		log_lines.append(entry["detail"])
		source_status.append(entry)
	elif bold_config is not None:
		# BOLDigger3 v3 (BOLD v5) classifies against the remote BOLD API and does
		# NOT use a local ``.ddb`` file, so we must NOT gate BOLD on ``.ddb``
		# presence — doing so silently skipped BOLD on every run. ``run_boldigger``
		# handles skipping cleanly (with an explicit reason surfaced to the user)
		# when the tool or credentials are missing.
		bold_result = run_boldigger_step(fasta_path, bold_config, output_dir / "bold")
		if bold_result.success:
			files.update(bold_result.files)
			bold_path = bold_result.files.get("bold_parquet") or bold_result.files.get("bold_xlsx")
			if isinstance(bold_path, Path) and bold_path.exists():
				# Top hit drives the consensus; the full hit list (when BOLD gave us
				# one) feeds the ranked table so the runner-up stays selectable.
				bold_all = bold_result.files.get("bold_all_hits_csv")
				_merge_rows(
					_collect_source_records("BOLD", bold_path),
					_collect_source_records("BOLD", bold_all) if isinstance(bold_all, Path) and bold_all.exists() else None,
				)
				# "ok" here means only that BOLD returned a file; whether it identified
				# anything is settled by _reconcile_source_status against the records.
				source_status.append({"source": "BOLD", "status": "ok", "detail": "identificação BOLD concluída", "key": "BOLD", "raw": bold_path})
			else:
				metrics["bold_error"] = 1.0
				log_lines.append("BOLD result missing parquet/xlsx path.")
				source_status.append({"source": "BOLD", "status": "failed", "detail": "resultado BOLD em falta"})
		else:
			metrics["bold_error"] = 1.0
			message = bold_result.error_message or "BOLD failed"
			log_lines.append(message)
			lowered = message.lower()
			status = "skipped" if ("credenciais" in lowered or "credentials" in lowered or "ignorado" in lowered) else "failed"
			source_status.append({"source": "BOLD", "status": status, "detail": message[:160]})

	# Persist a per-source status report so the UI can show, colour-coded, exactly
	# which reference databases ran, which found nothing, and which failed (e.g.
	# out-of-memory) — instead of a silent "None" that hides a crashed source.
	source_status_path = output_dir / "source_status.csv"
	if source_status:
		# Last word on every status: what the source actually put in the records, not
		# what its own proxy (a file existing, a BLAST hit count) suggested.
		source_status = _reconcile_source_status(source_status, records_by_asv, _count_fasta_records(fasta_path), marker=marker)
		_write_csv(source_status_path, source_status, ["source", "status", "detail"])
		files["source_status_csv"] = source_status_path

	if not records_by_asv:
		# A custom BLAST comparison is a valid standalone result even when no
		# consensus source (SILVA/PR2/NCBI/BOLD) produced records.
		if custom_blast_success:
			log_path.write_text("\n".join(log_lines) or "Custom BLAST comparison completed.", encoding="utf-8")
			return StepOutput(step_name="taxonomy", success=True, files=files, log=log_path, metrics=metrics)
		log_path.write_text("\n".join(log_lines) or "No taxonomy sources succeeded.", encoding="utf-8")
		return StepOutput(step_name="taxonomy", success=False, files=files, log=log_path, metrics=metrics, error_message="Nenhuma fonte de taxonomia teve sucesso.")

	consensus_rows = [_build_consensus_row(asv_id, records, consensus_method, feature_type) for asv_id, records in sorted(records_by_asv.items())]
	conflicts = _build_conflict_report(records_by_asv)
	conflict_columns = ["ASV_ID", "conflict_levels", "sources"]
	_write_csv(conflict_report, conflicts, conflict_columns)

	# Every candidate hit, ranked best-first per feature. This is the table the user
	# sorts, filters and picks the runner-up from.
	_write_csv(hits_path, _build_hits_rows(hits_by_asv, feature_type), HITS_COLUMNS)
	files["taxonomy_hits_csv"] = hits_path

	# One row per feature: the consensus call and the best hit of each evidence family,
	# in separate columns. Written here, from the same records, so the UI can render this
	# file rather than recompute it — one source, two consumers.
	summary_path = output_dir / "taxonomy_summary.csv"
	_write_csv(summary_path, _build_summary_rows(records_by_asv, consensus_rows, feature_type), SUMMARY_COLUMNS)
	files["taxonomy_summary_csv"] = summary_path

	# Per-source assignments carry the evidence (identity, coverage, confidence,
	# rank) that the collapsed consensus row throws away.
	assignments_path = output_dir / "source_assignments.csv"
	_write_csv(assignments_path, _build_assignment_rows(records_by_asv), ASSIGNMENT_COLUMNS)
	files["source_assignments_csv"] = assignments_path
	consensus_table = ConsensusTable(consensus_rows, list(CONSENSUS_COLUMNS))
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