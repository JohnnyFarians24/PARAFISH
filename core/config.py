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


def get_db_paths() -> DBPaths:
	"""Resolve database paths from environment variables."""
	return DBPaths(
		silva=_path_or_none(os.getenv("DB_SILVA_PATH", os.getenv("SILVA_FASTA_PATH", ""))),
		pr2=_path_or_none(os.getenv("DB_PR2_PATH", os.getenv("PR2_FASTA_PATH", ""))),
		pr2_dada2=_path_or_none(os.getenv("DB_PR2_DADA2_PATH", os.getenv("PR2_DADA2_FASTA_PATH", ""))),
		bold_dir=_path_or_none(os.getenv("DB_BOLD_PATH", os.getenv("BOLD_DB_PATH", ""))),
		blast=_path_or_none(os.getenv("DB_BLAST_PATH", os.getenv("NCBI_BLAST_DB", ""))),
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


def get_runs_dir() -> Path:
	"""Return the root directory for run inputs and outputs."""
	return Path(os.getenv("RUNS_DIR", "/data/runs").strip() or "/data/runs")


DB_ROOT = _resolve_db_root()
ALLOW_DB_AUTO_DOWNLOAD = _env_bool("ALLOW_DB_AUTO_DOWNLOAD", default=False)
KEEP_INTERMEDIATES = _env_bool("KEEP_INTERMEDIATES", default=False)
DELETE_AFTER_STEPS = [step.strip() for step in os.getenv("DELETE_AFTER_STEPS", "").split(",") if step.strip()]
BOLD_CACHE_ROOT = os.getenv("BOLD_CACHE_ROOT", os.path.join(DB_ROOT, "bold", "cache"))
os.makedirs(BOLD_CACHE_ROOT, exist_ok=True)
NCBI_REMOTE_FALLBACK = _env_bool("NCBI_REMOTE_FALLBACK", default=True)
BOLDIGGER_PYTHON = os.getenv("BOLDIGGER_PYTHON", os.getenv("BOLDIGGER_VENV_PYTHON", "")).strip()

BLAST_TIMEOUT_SEC = _env_int("BLAST_TIMEOUT_SEC", 900)
MAKEBLASTDB_TIMEOUT_SEC = _env_int("MAKEBLASTDB_TIMEOUT_SEC", 1800)
BOLDIGGER_TIMEOUT_SEC = _env_int("BOLDIGGER_TIMEOUT_SEC", 900)

# Timeout used specifically for remote NCBI BLAST calls. Remote NCBI BLAST
# is generally unsuitable for synchronous pipeline runs (it's slow and can
# hang the pipeline). Default to 0 (disabled) so users must opt-in or
# configure a local BLAST DB instead.
NCBI_TIMEOUT_SEC = _env_int("NCBI_TIMEOUT_SEC", 0)

FASTQC_PATH = _tool_cmd("fastqc", env_var="FASTQC_PATH")
VSEARCH_PATH = _tool_cmd("vsearch", env_var="VSEARCH_PATH")
MAFFT_PATH = _tool_cmd("mafft", "mafft.bat", env_var="MAFFT_PATH")
FASTTREE_PATH = _tool_cmd("FastTree", "fasttree", env_var="FASTTREE_PATH")
RAXML_PATH = _tool_cmd("raxml-ng", "raxmlng", env_var="RAXML_PATH")
PEAR_PATH = _tool_cmd("pear", env_var="PEAR_PATH")
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
