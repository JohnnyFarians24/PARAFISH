"""Central configuration constants and environment-driven runtime settings."""

from __future__ import annotations

import os
import shutil
import warnings as _warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


PIPELINE_STAGES = (
	"sample_selection",
	"qc",
	"preprocess",
	"analysis",
	"databases",
	"phylo",
)

SUPPORTED_MODES = ("SE", "PE")
SUPPORTED_MARKERS = ("COI", "18S", "ITS", "28S")

# Minimum allowed stage transitions for contiguous workflows.
ALLOWED_TRANSITIONS = {
	"sample_selection": ("qc",),
	"qc": ("preprocess",),
	"preprocess": ("analysis",),
	"analysis": ("databases", "phylo"),
	"databases": ("phylo",),
	"phylo": (),
}


@dataclass(frozen=True)
class DBPaths:
	silva: Optional[Path]
	pr2: Optional[Path]
	pr2_dada2: Optional[Path]
	bold_dir: Optional[Path]
	blast: Optional[Path]
	pr2_tax: Optional[Path] = None
	# SILVA's rank file (tax_slv_ssu_<release>.txt[.gz]). Opt-in: when unset, SILVA
	# contributes only domain and organism name, exactly as before.
	silva_ranks: Optional[Path] = None


@dataclass(frozen=True)
class FeatureFlags:
	enable_silva: bool
	enable_pr2: bool
	enable_bold: bool
	enable_ncbi: bool


@dataclass(frozen=True)
class BOLDCredentials:
	username: str
	password: str


def _env_bool(name, default=False):
	value = str(os.getenv(name, "")).strip().lower()
	if value in ("1", "true", "yes", "y", "on"):
		return True
	if value in ("0", "false", "no", "n", "off"):
		return False
	return bool(default)


def _env_int(name, default):
	raw = os.getenv(name, "").strip()
	if not raw:
		return int(default)
	try:
		value = int(raw)
	except Exception:
		return int(default)
	return value


def _load_env_file(path: Path) -> None:
	if not path.exists() or not path.is_file():
		return
	try:
		content = path.read_text(encoding="utf-8")
	except UnicodeError:
		content = path.read_text(encoding="utf-16")
	for raw_line in content.splitlines():
		line = raw_line.strip()
		if not line or line.startswith("#") or "=" not in line:
			continue
		key, value = line.split("=", 1)
		key = key.strip()
		value = value.strip().strip('"').strip("'")
		if key and key not in os.environ:
			os.environ[key] = value


def _bootstrap_env() -> None:
	_load_env_file(Path("/data/config/.env"))
	_load_env_file(Path.cwd() / ".env")


_bootstrap_env()


def _resolve_db_root():
	preferred = os.getenv("GLOBAL_DB_ROOT") or os.getenv("DB_ROOT") or "/data/dbs"
	fallback = os.path.join(os.getcwd(), "databases")
	try:
		os.makedirs(preferred, exist_ok=True)
		return preferred
	except Exception:
		os.makedirs(fallback, exist_ok=True)
		return fallback


def _tool_cmd(*names, env_var=None, fallback=None):
	if env_var:
		override = os.getenv(env_var)
		if override:
			return override
	for name in names:
		found = shutil.which(name)
		if found:
			return found
	if fallback:
		return fallback
	if names:
		return names[0]
	return ""


def _path_or_none(value: str) -> Optional[Path]:
	value = value.strip()
	if not value:
		return None
	return Path(value)


def _resolve_pr2_tax(pr2_fasta: Optional[Path]) -> Optional[Path]:
	"""Locate the PR2 mothur taxonomy file.

	Honour an explicit env override first; otherwise look next to the PR2 FASTA
	(the tax file normally ships alongside it) before falling back to the default
	database directory. This avoids silent "tax file not found" failures when the
	user points DB_PR2_PATH at a folder we don't otherwise know about.
	"""
	explicit = _path_or_none(os.getenv("DB_PR2_TAX_PATH", os.getenv("PR2_MOTHUR_TAX_PATH", "")))
	if explicit and explicit.exists():
		return explicit
	search_dirs = []
	if pr2_fasta is not None:
		search_dirs.append(pr2_fasta.parent)
	search_dirs.append(Path(PR2_DB_DIR))
	for directory in search_dirs:
		if not directory.exists():
			continue
		exact = directory / "pr2_version_5.1.1_SSU_mothur.tax.gz"
		if exact.exists():
			return exact
		for pattern in ("*mothur.tax.gz", "*mothur.tax", "*.tax.gz", "*.tax"):
			matches = sorted(directory.glob(pattern))
			if matches:
				return matches[0]
	return explicit or _path_or_none(DEFAULT_PR2_MOTHUR_TAX_PATH)


def get_db_paths() -> DBPaths:
	"""Resolve database paths from environment variables."""
	pr2 = _path_or_none(os.getenv("DB_PR2_PATH", os.getenv("PR2_FASTA_PATH", DEFAULT_PR2_MOTHUR_FASTA_PATH)))
	return DBPaths(
		silva=_path_or_none(os.getenv("DB_SILVA_PATH", os.getenv("SILVA_FASTA_PATH", ""))),
		pr2=pr2,
		pr2_dada2=_path_or_none(os.getenv("DB_PR2_DADA2_PATH", os.getenv("PR2_DADA2_FASTA_PATH", DEFAULT_PR2_DADA2_FASTA_PATH))),
		bold_dir=_path_or_none(os.getenv("DB_BOLD_PATH", os.getenv("BOLD_DB_PATH", ""))),
		blast=_path_or_none(os.getenv("DB_BLAST_PATH", os.getenv("NCBI_BLAST_DB", ""))),
		pr2_tax=_resolve_pr2_tax(pr2),
		silva_ranks=_path_or_none(os.getenv("DB_SILVA_RANKS_PATH", "")),
	)


def get_feature_flags() -> FeatureFlags:
	"""Return feature flags controlling optional sources and fallbacks."""
	return FeatureFlags(
		enable_silva=_env_bool("ENABLE_SILVA", default=True),
		enable_pr2=_env_bool("ENABLE_PR2", default=True),
		enable_bold=_env_bool("ENABLE_BOLD", default=True),
		enable_ncbi=_env_bool("ENABLE_NCBI", default=False),
	)


def get_bold_credentials() -> Optional[BOLDCredentials]:
	"""Read BOLD credentials from environment variables."""
	username = os.getenv("BOLD_USERNAME", "").strip()
	password = os.getenv("BOLD_PASSWORD", "").strip()
	if username and password:
		return BOLDCredentials(username=username, password=password)
	return None


def get_entrez_config() -> dict:
	"""NCBI Entrez identity for taxonomy lineage lookups.

	NCBI asks every automated caller to identify itself with an email; an optional
	API key raises the rate limit from 3 to 10 requests/second. Both are read from
	the environment so nothing is hard-coded, and both are optional — the lineage
	resolver degrades gracefully (upper ranks left blank, with a warning) when the
	email is missing or the network is down.
	"""
	return {
		"email": os.getenv("ENTREZ_EMAIL", "").strip(),
		"api_key": os.getenv("NCBI_API_KEY", os.getenv("ENTREZ_API_KEY", "")).strip(),
	}


def get_runs_dir() -> Path:
	"""Return the root directory for run inputs and outputs."""
	return Path(os.getenv("RUNS_DIR", "/data/runs").strip() or "/data/runs")


def check_core_tools() -> dict:
	"""Report availability of the external bioinformatics tools on PATH.

	The pipeline shells out to these tools directly (installed natively or inside
	Docker), so a live PATH check is a truthful health indicator — unlike a
	hard-coded "Docker offline" badge. Returns
	``{"tools": {name: bool, ...}, "available": int, "total": int, "ok": bool}``.
	"""
	checks = {
		"fastqc": FASTQC_PATH,
		"vsearch": VSEARCH_PATH,
		"blastn": BLASTN_PATH,
		"mafft": MAFFT_PATH,
		"fasttree": FASTTREE_PATH,
	}
	tools = {}
	for name, path in checks.items():
		tools[name] = bool(path and (shutil.which(path) or os.path.isfile(path)))
	available = sum(1 for present in tools.values() if present)
	total = len(tools)
	return {"tools": tools, "available": available, "total": total, "ok": available == total}


DB_ROOT = _resolve_db_root()
ALLOW_DB_AUTO_DOWNLOAD = _env_bool("ALLOW_DB_AUTO_DOWNLOAD", default=False)
KEEP_INTERMEDIATES = _env_bool("KEEP_INTERMEDIATES", default=False)
DELETE_AFTER_STEPS = [step.strip() for step in os.getenv("DELETE_AFTER_STEPS", "").split(",") if step.strip()]
BOLD_CACHE_ROOT = os.getenv("BOLD_CACHE_ROOT", os.path.join(DB_ROOT, "bold", "cache"))
os.makedirs(BOLD_CACHE_ROOT, exist_ok=True)
NCBI_REMOTE_FALLBACK = _env_bool("NCBI_REMOTE_FALLBACK", default=True)
BOLDIGGER_PYTHON = os.getenv("BOLDIGGER_PYTHON", os.getenv("BOLDIGGER_VENV_PYTHON", "")).strip()

BLAST_TIMEOUT_SEC = _env_int("BLAST_TIMEOUT_SEC", 900)

# vsearch --minuniquesize: the smallest cluster allowed to SEED an OTU during
# dereplication. 1 keeps singletons as possible seeds; 2 requires a sequence to be seen
# at least twice before it may found an OTU.
#
# Default 1, because this is a parasite-detection platform. The reads a floor of 2
# uniquely discards are not sequencing errors — those are within 3% of an abundant
# sequence and get absorbed into its OTU by clustering either way. What a floor of 2
# discards is the read that is globally unique AND >3% from everything else, which is the
# exact profile of a rare novel parasite. That false negative is the most expensive
# outcome this pipeline has, and it is not recoverable after the fact.
#
# Measured on the reference run (18S, primers removed): 1 and 2 produce IDENTICAL output
# — same 9 OTUs, same 154 223 reads mapped, zero 1-read OTUs, zero unmapped. The feared
# explosion of error-derived singleton OTUs does not happen. The only real cost is
# clustering time (57 s -> 247 s on that run), which scales with the number of uniques.
#
# So 2 remains a legitimate choice for a large host-dominated study buying back time, and
# is exposed rather than hardcoded: `read_tracking.csv` now carries `unmapped_reads`, so
# whatever a floor of 2 costs shows up as a number instead of disappearing. That is what
# makes it the user's call and not a code edit.
MIN_UNIQUE_SIZE = _env_int("PARAFISH_MIN_UNIQUE_SIZE", 1)

# Drop reads in which cutadapt found no primer, instead of passing them through uncut.
#
# On by default, and the reason is provenance rather than signal. A read with no primer
# match cannot be shown to be an amplicon of the targeted region — it may be off-target
# product, adapter-dimer, or a read whose 5' end was damaged. Keeping it means a feature
# in the results table that no one can attest came from the assay that was run.
#
# Measured on the COI/Leray run (aptl64b_01, 207 528 reads): 92.32% carried the forward
# primer, so cutadapt discards ~7.7% of INPUT reads here. The effect on the RESULT is
# far smaller — only 335 of 179 054 filtered reads (0.19%) currently survive to
# clustering at an untrimmed length, because PEAR assembly and quality filtering already
# remove most of the rest. Expect read_tracking's input->trimmed drop to widen by ~7.7
# points; that is this flag working, not a regression.
#
# Independent of MIN_UNIQUE_SIZE, and not in tension with it: that one keeps genuinely
# rare amplicons, this one removes material that is not a verified amplicon at all.
DISCARD_UNTRIMMED = _env_bool("PARAFISH_DISCARD_UNTRIMMED", True)

# Use a LOCAL BLAST database instead of NCBI's public remote service. Off by default:
# the public remote nt is the source the project's supervisors asked for, and it is what
# a reviewer can reproduce without our filesystem. A local db is the emergency exit for
# when NCBI is down or refusing — and, being our own CPU, the one place a sensitive
# search (task=blastn, word_size 11) is affordable; the public quota refuses that
# outright (SIGXCPU). Set PARAFISH_NCBI_LOCAL_DB=1 and point BLASTDB/blast db paths at it.
NCBI_USE_LOCAL_DB = _env_bool("PARAFISH_NCBI_LOCAL_DB", default=False)
MAKEBLASTDB_TIMEOUT_SEC = _env_int("MAKEBLASTDB_TIMEOUT_SEC", 1800)
BOLDIGGER_TIMEOUT_SEC = _env_int("BOLDIGGER_TIMEOUT_SEC", 900)

# Timeout for remote NCBI BLAST. Remote BLAST runs on NCBI's servers, so it costs
# no local RAM, but it is slow (~1 min for a single query, and NCBI queues under
# load). 1800s is enough for a typical OTU set; raise it for large ones.
NCBI_TIMEOUT_SEC = _env_int("NCBI_TIMEOUT_SEC", 1800)

# Directory holding NCBI's taxdb (taxdb.btd/.bti). Without it, blastn cannot turn a
# taxid into a species name and every `sscinames` column comes back "N/A" — the hit
# is there but unusable for taxonomy. Exported into the blastn environment.
BLASTDB_PATH = os.getenv("BLASTDB", "").strip()

FASTQC_PATH = _tool_cmd("fastqc", env_var="FASTQC_PATH")
VSEARCH_PATH = _tool_cmd("vsearch", env_var="VSEARCH_PATH")
MAFFT_PATH = _tool_cmd("mafft", "mafft.bat", env_var="MAFFT_PATH")
FASTTREE_PATH = _tool_cmd("FastTree", "fasttree", env_var="FASTTREE_PATH")
RAXML_PATH = _tool_cmd("raxml-ng", "raxmlng", env_var="RAXML_PATH")
PEAR_PATH = _tool_cmd("pear", env_var="PEAR_PATH")
CUTADAPT_PATH = _tool_cmd("cutadapt", env_var="CUTADAPT_PATH")
BLASTN_PATH = _tool_cmd("blastn", env_var="BLASTN_PATH")
MAKEBLASTDB_PATH = _tool_cmd("makeblastdb", env_var="MAKEBLASTDB_PATH")
RSCRIPT_PATH = _tool_cmd("Rscript", env_var="RSCRIPT_PATH")

TRIMMOMATIC_JAR_PATH = os.getenv(
	"TRIMMOMATIC_JAR_PATH",
	"/opt/Trimmomatic/trimmomatic.jar",
)
TRIMMOMATIC_ADAPTERS = os.getenv(
	"TRIMMOMATIC_ADAPTERS",
	"/opt/Trimmomatic/adapters/TruSeq3-SE.fa",
)
TRIMMOMATIC_ADAPTERS_PE = os.getenv(
	"TRIMMOMATIC_ADAPTERS_PE",
	"/opt/Trimmomatic/adapters/TruSeq3-PE.fa",
)

if not os.path.exists(TRIMMOMATIC_JAR_PATH):
	_warnings.warn(
		f"Trimmomatic JAR not found: {TRIMMOMATIC_JAR_PATH}. "
		"Define TRIMMOMATIC_JAR_PATH in the environment.",
		RuntimeWarning,
		stacklevel=1,
	)

DEFAULT_SILVA_URL = os.getenv(
	"SILVA_FASTA_URL",
	"https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz",
)
DEFAULT_SILVA_PATH = os.getenv("SILVA_FASTA_PATH", "")
DEFAULT_PR2_URL = os.getenv(
	"PR2_FASTA_URL",
	"https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_mothur.fasta.gz",
)
DEFAULT_PR2_PATH = os.getenv("PR2_FASTA_PATH", "")
PR2_DB_DIR = os.path.join(DB_ROOT, "pr2", "5_1_1")
DEFAULT_PR2_MOTHUR_FASTA_PATH = os.getenv(
	"PR2_MOTHUR_FASTA_PATH",
	os.path.join(PR2_DB_DIR, "pr2_version_5.1.1_SSU_mothur.fasta.gz"),
)
DEFAULT_PR2_MOTHUR_TAX_PATH = os.getenv(
	"PR2_MOTHUR_TAX_PATH",
	os.path.join(PR2_DB_DIR, "pr2_version_5.1.1_SSU_mothur.tax.gz"),
)
DEFAULT_PR2_DADA2_FASTA_PATH = os.getenv(
	"PR2_DADA2_FASTA_PATH",
	os.path.join(PR2_DB_DIR, "pr2_version_5.1.1_SSU_dada2.fasta.gz"),
)
NCBI_BLAST_DB = os.getenv("NCBI_BLAST_DB", "")
