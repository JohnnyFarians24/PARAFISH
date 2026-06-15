import base64
import dash
from dash import dcc, html, Input, Output, State, dash_table
import dash_bootstrap_components as dbc
from datetime import datetime
import importlib.util
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import zipfile
import traceback
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import dash_cytoscape as cyto
import dash_uploader as du
from dash.exceptions import PreventUpdate
from Bio import Phylo
from flask import send_from_directory
from scipy.stats import entropy
from scipy.spatial.distance import pdist, squareform
from sample_parsing import parse_fastq_name, validate_manifest_rows
import jobs
import vsearch_pipeline as vp

# --- Monkey-patch dash_uploader for Dash 3.x compatibility ---
import dash_uploader.callbacks as _du_callbacks
_original_create_dash_callback = _du_callbacks.create_dash_callback

def _fixed_create_dash_callback(callback, settings_ref):
    original_wrapper = _original_create_dash_callback(callback, settings_ref)
    def safe_wrapper(iscompleted, filenames, upload_id):
        if not iscompleted:
            raise PreventUpdate
        return original_wrapper(iscompleted, filenames, upload_id)
    return safe_wrapper

_du_callbacks.create_dash_callback = _fixed_create_dash_callback
# --- End monkey-patch ---

UPLOAD_DIRECTORY = os.path.join(os.getcwd(), "uploads")
os.makedirs(UPLOAD_DIRECTORY, exist_ok=True)

PROJECTS_ROOT = os.path.join(os.getcwd(), "projects")
os.makedirs(PROJECTS_ROOT, exist_ok=True)

_LAST_EXPORT_ERROR = ""


def _project_exports_dir(project_dir):
    if not project_dir:
        return None
    exports_dir = os.path.join(project_dir, "exports")
    os.makedirs(exports_dir, exist_ok=True)
    return exports_dir


def _set_export_error(message=""):
    global _LAST_EXPORT_ERROR
    _LAST_EXPORT_ERROR = str(message or "")


def _is_export_blocked_component(component):
    lowered = str(component or "").lower()
    return lowered in {"temp", "tmp", "intermediate"} or lowered.startswith("temp_") or lowered.startswith("tmp_")


def _should_include_export_file(file_path, run_dir):
    if not file_path or not run_dir or not os.path.isfile(file_path):
        return False
    if os.path.getsize(file_path) <= 0:
        return False

    base_name = os.path.basename(file_path).lower()
    if base_name.endswith((".zip", ".tmp", ".temp")):
        return False
    if base_name.startswith(("temp_", "tmp_")):
        return False

    rel_path = os.path.relpath(file_path, run_dir)
    if rel_path.startswith(".."):
        return False

    parts = rel_path.split(os.sep)
    if not parts or parts[0] not in {"steps", "results"}:
        return False
    if any(_is_export_blocked_component(part) for part in parts):
        return False
    return True


def _env_bool(name, default=False):
    value = str(os.getenv(name, "")).strip().lower()
    if value in ("1", "true", "yes", "y", "on"):
        return True
    if value in ("0", "false", "no", "n", "off"):
        return False
    return bool(default)


def _resolve_db_root():
    preferred = os.getenv("GLOBAL_DB_ROOT") or os.getenv("DB_ROOT") or "/data/dbs"
    fallback = os.path.join(os.getcwd(), "databases")
    try:
        os.makedirs(preferred, exist_ok=True)
        return preferred
    except Exception:
        os.makedirs(fallback, exist_ok=True)
        return fallback


DB_ROOT = _resolve_db_root()
ALLOW_DB_AUTO_DOWNLOAD = _env_bool("ALLOW_DB_AUTO_DOWNLOAD", default=False)
KEEP_INTERMEDIATES = _env_bool("KEEP_INTERMEDIATES", default=False)
DELETE_AFTER_STEPS = [step.strip() for step in os.getenv("DELETE_AFTER_STEPS", "").split(",") if step.strip()]
BOLD_CACHE_ROOT = os.getenv("BOLD_CACHE_ROOT", os.path.join(DB_ROOT, "bold", "cache"))
os.makedirs(BOLD_CACHE_ROOT, exist_ok=True)
NCBI_REMOTE_FALLBACK = _env_bool("NCBI_REMOTE_FALLBACK", default=True)
BOLDIGGER_PYTHON = os.getenv("BOLDIGGER_PYTHON", os.getenv("BOLDIGGER_VENV_PYTHON", "")).strip()


def _env_int(name, default):
    raw = os.getenv(name, "").strip()
    if not raw:
        return int(default)
    try:
        value = int(raw)
    except Exception:
        return int(default)
    return value


BLAST_TIMEOUT_SEC = _env_int("BLAST_TIMEOUT_SEC", 900)
MAKEBLASTDB_TIMEOUT_SEC = _env_int("MAKEBLASTDB_TIMEOUT_SEC", 1800)
BOLDIGGER_TIMEOUT_SEC = _env_int("BOLDIGGER_TIMEOUT_SEC", 900)

TRIMMOMATIC_JAR = os.path.join(os.getcwd(), "Trimmomatic-0.39", "trimmomatic-0.39.jar")
TRIMMOMATIC_ADAPTERS = os.path.join(os.getcwd(), "Trimmomatic-0.39", "adapters", "TruSeq3-SE.fa")
TRIMMOMATIC_ADAPTERS_PE = os.path.join(os.getcwd(), "Trimmomatic-0.39", "adapters", "TruSeq3-PE.fa")


app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True
)
server = app.server
server.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024 * 1024
du.configure_upload(app, UPLOAD_DIRECTORY)


@server.route('/jobs/<run_id>/status')
def _job_status_endpoint(run_id):
    try:
        meta = jobs.job_status(run_id)
        if meta is None:
            return ("Not Found", 404)
        return json.dumps(meta), 200, {"Content-Type": "application/json"}
    except Exception:
        return ("Internal Error", 500)


@server.route('/jobs/<run_id>/log')
def _job_log_endpoint(run_id):
    try:
        log = jobs.job_log(run_id)
        if log is None:
            return ("Not Found", 404)
        return log, 200, {"Content-Type": "text/plain"}
    except Exception:
        return ("Internal Error", 500)


# ============================================================================
# INITIALIZE STORAGE CLEANUP
# ============================================================================
def _get_all_project_dirs():
    """Get list of all project directories."""
    projects = []
    if not os.path.exists(PROJECTS_ROOT):
        return projects
    try:
        for entry in os.listdir(PROJECTS_ROOT):
            pdir = os.path.join(PROJECTS_ROOT, entry)
            if os.path.isdir(pdir):
                projects.append(pdir)
    except Exception:
        pass
    return projects

def _tool_cmd(name, env_var=None, fallback=None):
    if env_var:
        override = os.getenv(env_var)
        if override:
            return override
    found = shutil.which(name)
    if found:
        return found
    return fallback or name


FASTQC_PATH = _tool_cmd("fastqc", "FASTQC_PATH")
VSEARCH_PATH = _tool_cmd("vsearch", "VSEARCH_PATH")
MAFFT_PATH = _tool_cmd("mafft", "MAFFT_PATH")
FASTTREE_PATH = _tool_cmd("fasttree", "FASTTREE_PATH")
RAXML_PATH = _tool_cmd("raxml-ng", "RAXML_PATH")
PEAR_PATH = _tool_cmd("pear", "PEAR_PATH")
BLASTN_PATH = _tool_cmd("blastn", "BLASTN_PATH")
MAKEBLASTDB_PATH = _tool_cmd("makeblastdb", "MAKEBLASTDB_PATH")
RSCRIPT_PATH = _tool_cmd("Rscript", "RSCRIPT_PATH")

DEFAULT_SILVA_URL = os.getenv(
    "SILVA_FASTA_URL",
    "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
)
DEFAULT_SILVA_PATH = os.getenv("SILVA_FASTA_PATH", "")
DEFAULT_PR2_URL = os.getenv(
    "PR2_FASTA_URL",
    "https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_mothur.fasta.gz"
)
DEFAULT_PR2_PATH = os.getenv("PR2_FASTA_PATH", "")
PR2_DB_DIR = os.path.join(DB_ROOT, "pr2", "5_1_1")
DEFAULT_PR2_MOTHUR_FASTA_PATH = os.getenv(
    "PR2_MOTHUR_FASTA_PATH",
    os.path.join(PR2_DB_DIR, "pr2_version_5.1.1_SSU_mothur.fasta.gz")
)
DEFAULT_PR2_MOTHUR_TAX_PATH = os.getenv(
    "PR2_MOTHUR_TAX_PATH",
    os.path.join(PR2_DB_DIR, "pr2_version_5.1.1_SSU_mothur.tax.gz")
)
DEFAULT_PR2_DADA2_FASTA_PATH = os.getenv(
    "PR2_DADA2_FASTA_PATH",
    os.path.join(PR2_DB_DIR, "pr2_version_5.1.1_SSU_dada2.fasta.gz")
)
NCBI_BLAST_DB = os.getenv("NCBI_BLAST_DB", "")


def _safe_name(text):
    value = str(text or "").strip()
    if not value:
        return "Project"
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    cleaned = cleaned.strip("._")
    return cleaned or "Project"


def _project_dir(project_name):
    return os.path.join(PROJECTS_ROOT, _safe_name(project_name))


def _ensure_project(project_name):
    safe_project = _safe_name(project_name)
    pdir = _project_dir(safe_project)
    os.makedirs(pdir, exist_ok=True)
    os.makedirs(os.path.join(pdir, "runs"), exist_ok=True)
    os.makedirs(os.path.join(pdir, "uploads"), exist_ok=True)
    os.makedirs(_project_exports_dir(pdir), exist_ok=True)
    meta_path = os.path.join(pdir, "project.json")
    if not os.path.exists(meta_path):
        _write_json(meta_path, {
            "project_name": safe_project,
            "created_at": datetime.now().isoformat(timespec="seconds")
        })
    return pdir


def _create_run(project_dir):
    runs_dir = os.path.join(project_dir, "runs")
    os.makedirs(runs_dir, exist_ok=True)
    stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    run_dir = os.path.join(runs_dir, stamp)
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


# ============================================================================
# STORAGE MANAGEMENT SYSTEM - Prevent infinite disk growth
# ============================================================================

class StorageConfig:
    """Configuration for storage management policies."""
    # Maximum number of old runs to keep per project
    MAX_RUNS_TO_KEEP = int(os.getenv("MAX_RUNS_TO_KEEP", "10"))
    # Maximum age of runs in days before cleanup
    MAX_RUN_AGE_DAYS = int(os.getenv("MAX_RUN_AGE_DAYS", "30"))
    # Enable automatic cleanup of temp files
    AUTO_DELETE_TEMP = _env_bool("AUTO_DELETE_TEMP", default=True)
    # Maximum temp file age in days
    TEMP_FILE_MAX_AGE_DAYS = int(os.getenv("TEMP_FILE_MAX_AGE_DAYS", "3"))
    # Delete these steps' intermediates after completion
    DELETE_INTERMEDIATE_STEPS = [s.strip() for s in os.getenv("DELETE_INTERMEDIATE_STEPS", "qc,preprocess").split(",") if s.strip()]


class StorageManager:
    """Manages storage cleanup and retention policies."""
    
    @staticmethod
    def cleanup_old_runs(project_dir):
        """Delete old runs, keeping only MAX_RUNS_TO_KEEP most recent."""
        runs_dir = os.path.join(project_dir, "runs")
        if not os.path.exists(runs_dir):
            return []
        
        deleted = []
        try:
            entries = []
            for entry in os.listdir(runs_dir):
                entry_path = os.path.join(runs_dir, entry)
                if os.path.isdir(entry_path):
                    mtime = os.path.getmtime(entry_path)
                    entries.append((mtime, entry_path, entry))
            
            entries.sort(reverse=True)
            
            # Keep only MAX_RUNS_TO_KEEP
            for mtime, entry_path, entry in entries[StorageConfig.MAX_RUNS_TO_KEEP:]:
                try:
                    shutil.rmtree(entry_path)
                    deleted.append(entry)
                except Exception as e:
                    pass
        except Exception:
            pass
        
        return deleted
    
    @staticmethod
    def cleanup_temp_files(project_dir):
        """Delete temporary files older than TEMP_FILE_MAX_AGE_DAYS."""
        deleted = []
        temp_patterns = [
            ("steps/*/temp_*", "step temp files"),
            ("runs/*/steps/*/tmp_*", "run temp files"),
            ("uploads/*_temp_*", "upload temp files"),
        ]
        
        now = time.time()
        max_age_seconds = StorageConfig.TEMP_FILE_MAX_AGE_DAYS * 86400
        
        try:
            import glob
            for pattern, desc in temp_patterns:
                full_pattern = os.path.join(project_dir, pattern)
                for filepath in glob.glob(full_pattern):
                    try:
                        file_age = now - os.path.getmtime(filepath)
                        if file_age > max_age_seconds:
                            if os.path.isfile(filepath):
                                os.remove(filepath)
                            elif os.path.isdir(filepath):
                                shutil.rmtree(filepath)
                            deleted.append(os.path.basename(filepath))
                    except Exception:
                        pass
        except Exception:
            pass
        
        return deleted
    
    @staticmethod
    def compact_run_outputs(run_dir):
        """Compress large output files to save space."""
        compressed_files = []
        try:
            large_patterns = [
                "steps/*/intermediate/*.fastq",
                "steps/*/intermediate/*.fq",
            ]
            import glob
            for pattern in large_patterns:
                full_pattern = os.path.join(run_dir, pattern)
                for filepath in glob.glob(full_pattern):
                    try:
                        if os.path.getsize(filepath) > 10 * 1024 * 1024:  # > 10MB
                            gz_path = filepath + ".gz"
                            if not os.path.exists(gz_path):
                                with open(filepath, 'rb') as f_in:
                                    with gzip.open(gz_path, 'wb') as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                os.remove(filepath)
                                compressed_files.append(os.path.basename(filepath))
                    except Exception:
                        pass
        except Exception:
            pass
        
        return compressed_files

    @staticmethod
    def cleanup_export_archives(project_dir, run_dir=None):
        """Remove stale export archives and stray run-level zip files."""
        deleted = []
        try:
            if run_dir and os.path.isdir(run_dir):
                export_staging_dir = os.path.join(run_dir, "export_staging")
                if os.path.isdir(export_staging_dir):
                    shutil.rmtree(export_staging_dir, ignore_errors=True)

                for name in os.listdir(run_dir):
                    file_path = os.path.join(run_dir, name)
                    if os.path.isfile(file_path) and name.lower().endswith(".zip"):
                        try:
                            os.remove(file_path)
                            deleted.append(file_path)
                        except OSError:
                            pass

            exports_dir = _project_exports_dir(project_dir) if project_dir else None
            if exports_dir and os.path.isdir(exports_dir):
                zip_entries = []
                for name in os.listdir(exports_dir):
                    file_path = os.path.join(exports_dir, name)
                    if os.path.isfile(file_path) and name.lower().endswith(".zip"):
                        try:
                            zip_entries.append((os.path.getmtime(file_path), file_path))
                        except OSError:
                            continue

                zip_entries.sort(key=lambda item: item[0], reverse=True)
                for _, file_path in zip_entries[StorageConfig.MAX_EXPORTS_TO_KEEP:]:
                    try:
                        os.remove(file_path)
                        deleted.append(file_path)
                    except OSError:
                        pass
        except Exception:
            pass

        return deleted


# Call cleanup on startup
try:
    for project_dir in _get_all_project_dirs():
        try:
            if StorageConfig.AUTO_DELETE_TEMP:
                StorageManager.cleanup_temp_files(project_dir)
            StorageManager.cleanup_old_runs(project_dir)
            StorageManager.cleanup_export_archives(project_dir)
        except Exception:
            pass
except Exception:
    pass


# ============================================================================
# PHASE 2: RUN EXECUTION INTEGRATION - Cleanup, PE handling, Export
# ============================================================================

def _apply_cleanup_post_run(project_dir, run_dir):
    """Apply cleanup after a run completes. Wraps StorageManager operations."""
    try:
        # Cleanup old runs (keep most recent)
        StorageManager.cleanup_old_runs(project_dir)
        
        # Cleanup temp files
        if StorageConfig.AUTO_DELETE_TEMP:
            StorageManager.cleanup_temp_files(project_dir)
        
        # Compact outputs in this run
        StorageManager.compact_run_outputs(run_dir)

        # Remove old export archives and stray run-level zip files
        StorageManager.cleanup_export_archives(project_dir, run_dir)
        
        return True
    except Exception as e:
        # Log but don't crash - cleanup is best-effort
        return False


def _ensure_pe_sample_paths(sample_info):
    """Ensure sample has correct r1/r2 or se paths for PE/SE distinction.
    
    Returns: Updated sample_info dict with validated paths
    """
    if not sample_info or not isinstance(sample_info, dict):
        return sample_info or {}
    
    info = sample_info.copy()
    paths = info.get('paths', {})
    kind = info.get('kind', 'SE')
    
    # If kind is PE, ensure both r1 and r2 exist
    if kind == 'PE':
        if not (paths.get('r1') and paths.get('r2')):
            # Downgrade to SE if r1/r2 not both present
            info['kind'] = 'SE'
            if paths.get('se'):
                pass  # Keep SE path
            elif paths.get('r1'):
                info['paths'] = {'se': paths['r1']}
            elif paths.get('r2'):
                info['paths'] = {'se': paths['r2']}
    
    return info


def _generate_run_export_zip(project_state, run_config, samples_store, workflow_scope):
    """Compatibility wrapper for the safe run export builder."""
    return _create_session_results_zip(project_state)


def _get_sample_fastq_pair(sample_id, samples_store):
    """Get r1/r2 file paths for a PE sample, or se path for SE sample.
    
    Returns: {'kind': 'PE'|'SE', 'files': [path(s)]}
    """
    samples = samples_store or {}
    sample = samples.get(sample_id, {})
    
    info = _ensure_pe_sample_paths(sample)
    paths = info.get('paths', {})
    kind = info.get('kind', 'SE')
    
    if kind == 'PE' and paths.get('r1') and paths.get('r2'):
        return {'kind': 'PE', 'files': [paths['r1'], paths['r2']]}
    elif paths.get('se'):
        return {'kind': 'SE', 'files': [paths['se']]}
    elif paths.get('r1'):
        return {'kind': 'SE', 'files': [paths['r1']]}
    elif paths.get('r2'):
        return {'kind': 'SE', 'files': [paths['r2']]}
    
    return {'kind': 'SE', 'files': []}


# ============================================================================
# SAMPLE GROUPING SYSTEM - Normalize and group PE pairs
# ============================================================================

def _detect_pe_pair_suffix(filename):
    """Detect if filename has PE pair suffix (_1/_2, R1/R2, .R1/.R2, etc.)."""
    lower = filename.lower()
    
    # Common PE patterns
    patterns = [
        (r'_1(?:\.\w+)?$', '1'),  # _1
        (r'_2(?:\.\w+)?$', '2'),  # _2
        (r'\.r1(?:\.\w+)?$', 'r1'),  # .R1
        (r'\.r2(?:\.\w+)?$', 'r2'),  # .R2
        (r'_r1(?:\.\w+)?$', 'r1'),  # _R1
        (r'_r2(?:\.\w+)?$', 'r2'),  # _R2
        (r'\.raw_1(?:\.\w+)?$', '1'),  # .raw_1
        (r'\.raw_2(?:\.\w+)?$', '2'),  # .raw_2
        (r'_raw_1(?:\.\w+)?$', '1'),  # _raw_1
        (r'_raw_2(?:\.\w+)?$', '2'),  # _raw_2
    ]
    
    for pattern, suffix_type in patterns:
        if re.search(pattern, lower):
            return suffix_type
    
    return None


def _normalize_sample_name(filename):
    """Extract base sample name from filename, removing PE suffixes and extensions."""
    name = os.path.basename(str(filename or "")).strip()
    if not name:
        return ""

    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fasta.gz", ".fa.gz", ".fasta", ".fa", ".gz"):
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    else:
        name = os.path.splitext(name)[0]

    suffix_patterns = [
        r"(?i)(?:[_\.-]?(?:raw)?[_\.-]?(?:r|fwd|rev)?[_\.-]?(?:1|2))$",
        r"(?i)(?:[_\.-](?:r1|r2))$",
        r"(?i)(?:[_\.-](?:forward|reverse))$",
    ]

    changed = True
    while changed and name:
        changed = False
        for pattern in suffix_patterns:
            new_name = re.sub(pattern, "", name)
            if new_name != name:
                name = new_name
                changed = True

    return name.strip('_.-')


def _is_paired_end_pair(file1, file2):
    """Check if two filenames represent a PE pair."""
    if not file1 or not file2:
        return False
    
    norm1 = _normalize_sample_name(file1)
    norm2 = _normalize_sample_name(file2)
    
    # Same base name = likely a pair
    return norm1 == norm2


def _group_pe_pairs(file_list):
    """Group files into (SE) and paired (R1, R2) groups. Returns list of (type, [files])."""
    if not file_list:
        return []
    
    processed = set()
    groups = []
    
    # First pass: detect pairs
    for i, file1 in enumerate(file_list):
        if i in processed:
            continue
        
        suffix1 = _detect_pe_pair_suffix(file1)
        if not suffix1:
            # Likely SE
            groups.append(("SE", [file1]))
            processed.add(i)
            continue
        
        # Look for pair
        found_pair = False
        for j, file2 in enumerate(file_list[i+1:], start=i+1):
            if j in processed:
                continue
            
            suffix2 = _detect_pe_pair_suffix(file2)
            if not suffix2:
                continue
            
            # Check if suffixes complement each other (1 with 2, R1 with R2)
            if (suffix1 in ('1', 'r1') and suffix2 in ('2', 'r2')) or \
               (suffix1 in ('2', 'r2') and suffix2 in ('1', 'r1')):
                if _is_paired_end_pair(file1, file2):
                    # Found a pair!
                    r1 = file1 if suffix1 in ('1', 'r1') else file2
                    r2 = file2 if suffix1 in ('1', 'r1') else file1
                    groups.append(("PE", [r1, r2]))
                    processed.add(i)
                    processed.add(j)
                    found_pair = True
                    break
        
        if not found_pair:
            # No pair found, treat as SE
            groups.append(("SE", [file1]))
            processed.add(i)
    
    return groups


def _get_active_root(project_state):
    state = project_state or {}
    run_dir = state.get("run_dir")
    if run_dir:
        os.makedirs(run_dir, exist_ok=True)
        return run_dir
    project_dir = state.get("project_dir")
    if project_dir:
        os.makedirs(project_dir, exist_ok=True)
        return project_dir
    return UPLOAD_DIRECTORY


def _is_within_directory(root_dir, target_path):
    root = os.path.abspath(root_dir)
    target = os.path.abspath(target_path)
    return target == root or target.startswith(root + os.sep)


def _safe_extract_zip(zip_path, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        for member in zf.infolist():
            member_name = member.filename
            if not member_name or member_name.endswith("/"):
                continue
            final_path = os.path.abspath(os.path.join(dest_dir, member_name))
            if not _is_within_directory(dest_dir, final_path):
                raise ValueError(f"Unsafe ZIP path detected: {member_name}")
        zf.extractall(dest_dir)


def _safe_extract_archive(archive_path, dest_dir):
    lower = archive_path.lower()
    if lower.endswith(".zip"):
        _safe_extract_zip(archive_path, dest_dir)
        return
    if lower.endswith(".rar"):
        os.makedirs(dest_dir, exist_ok=True)
        unrar_bin = shutil.which("unrar")
        if unrar_bin:
            cmd = [unrar_bin, "x", "-y", archive_path, dest_dir]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode == 0:
                return
        bsdtar_bin = shutil.which("bsdtar")
        if bsdtar_bin:
            cmd = [bsdtar_bin, "-xf", archive_path, "-C", dest_dir]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode == 0:
                return
        raise RuntimeError("RAR extraction failed. Install 'unrar' or 'bsdtar'.")
    raise ValueError("Unsupported archive format. Please upload ZIP or RAR.")


def _index_samples_from_dir(dataset_dir):
    candidates = []
    for root, _, files in os.walk(dataset_dir):
        for name in files:
            lower = name.lower()
            if lower.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fasta", ".fa")):
                candidates.append(os.path.join(root, name))
    return assign_samples_from_files(sorted(candidates), current_samples={})

def _download_href(file_path):
    file_path = os.path.abspath(file_path)
    if file_path.startswith(os.path.abspath(UPLOAD_DIRECTORY) + os.sep):
        rel = os.path.relpath(file_path, UPLOAD_DIRECTORY).replace(os.sep, "/")
        return f"/uploads/{rel}"
    cwd = os.path.abspath(os.getcwd()) + os.sep
    if file_path.startswith(cwd):
        rel = os.path.relpath(file_path, os.getcwd()).replace(os.sep, "/")
        return f"/files/{rel}"
    return "#"

def _step_base(project_state, step_name):
    base = _get_active_root(project_state)
    if base == UPLOAD_DIRECTORY:
        return base
    step_dir = os.path.join(base, "steps", step_name)
    os.makedirs(step_dir, exist_ok=True)
    return step_dir


def _step_outputs_base(project_state, step_name):
    step_dir = _step_base(project_state, step_name)
    outputs_dir = os.path.join(step_dir, "outputs")
    os.makedirs(outputs_dir, exist_ok=True)
    return outputs_dir


def _step_output_path(project_state, step_name, filename, sample_id=None):
    outputs_dir = _step_outputs_base(project_state, step_name)
    if sample_id:
        outputs_dir = os.path.join(outputs_dir, _safe_name(sample_id))
        os.makedirs(outputs_dir, exist_ok=True)
    return os.path.join(outputs_dir, filename) if filename else outputs_dir


def _mirror_output_to_results(source_path, results_path):
    if not source_path or not results_path:
        return results_path
    if not os.path.exists(source_path) or os.path.getsize(source_path) == 0:
        return results_path
    os.makedirs(os.path.dirname(results_path), exist_ok=True)
    try:
        if os.path.abspath(source_path) != os.path.abspath(results_path):
            shutil.copy2(source_path, results_path)
    except Exception:
        return results_path
    return results_path


def _write_text_file(path, content):
    if not path:
        return None
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(content if content is not None else "")
    return path


def _discover_dada2_preprocess_pairs(project_state):
    preprocess_root = _step_base(project_state, "preprocess")
    if not os.path.exists(preprocess_root):
        return []

    pairs = {}
    for root, _, files in os.walk(preprocess_root):
        for name in files:
            lower = name.lower()
            path = os.path.join(root, name)
            if not lower.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
                continue

            if "_trimmed_r1_paired" in lower:
                key = re.sub(r"(?i)_trimmed_r1_paired(?:\.fastq(?:\.gz)?|\.fq(?:\.gz)?)$", "", name)
                pairs.setdefault(key, {})["r1"] = path
            elif "_trimmed_r2_paired" in lower:
                key = re.sub(r"(?i)_trimmed_r2_paired(?:\.fastq(?:\.gz)?|\.fq(?:\.gz)?)$", "", name)
                pairs.setdefault(key, {})["r2"] = path

    resolved = []
    for key, item in sorted(pairs.items()):
        r1 = item.get("r1")
        r2 = item.get("r2")
        if r1 and r2 and _path_is_valid(r1) and _path_is_valid(r2):
            resolved.append((key, r1, r2))
    return resolved


def _resolve_dada2_inputs(project_state, uploaded_files, samples_store, workflow_scope):
    selected_samples = list(_iter_workflow_samples(samples_store, workflow_scope))
    if not selected_samples:
        return {"mode": None, "fwd_files": [], "rev_files": [], "single_files": [], "source": None, "error": "No samples selected for DADA2."}

    pe_entries = []
    se_entries = []
    for sid, sample in selected_samples:
        info = sample or {}
        paths = info.get("paths", {}) if isinstance(info, dict) else {}
        sample_mode = str(info.get("mode") or info.get("kind") or "SE").upper()

        if sample_mode == "PE":
            r1 = paths.get("r1")
            r2 = paths.get("r2")
            if not (r1 and r2 and _path_is_valid(r1) and _path_is_valid(r2)):
                return {
                    "mode": None,
                    "fwd_files": [],
                    "rev_files": [],
                    "single_files": [],
                    "source": "workflow_scope",
                    "error": f"PE sample {sid} is missing a valid forward/reverse pair.",
                }
            pe_entries.append((sid, r1, r2))
        else:
            se = paths.get("se")
            if not (se and _path_is_valid(se)):
                return {
                    "mode": None,
                    "fwd_files": [],
                    "rev_files": [],
                    "single_files": [],
                    "source": "workflow_scope",
                    "error": f"SE sample {sid} is missing a valid single-end file.",
                }
            se_entries.append((sid, se))

    if pe_entries and se_entries:
        return {
            "mode": None,
            "fwd_files": [],
            "rev_files": [],
            "single_files": [],
            "source": "workflow_scope",
            "error": "Mixed SE and PE samples are not supported in a single DADA2 run.",
        }

    if pe_entries:
        pe_entries = sorted(pe_entries, key=lambda item: item[0])
        return {
            "mode": "PE",
            "fwd_files": [item[1] for item in pe_entries],
            "rev_files": [item[2] for item in pe_entries],
            "single_files": [],
            "source": "workflow_scope_pairs",
            "error": "",
        }

    if se_entries:
        se_entries = sorted(se_entries, key=lambda item: item[0])
        return {
            "mode": "SE",
            "fwd_files": [],
            "rev_files": [],
            "single_files": [item[1] for item in se_entries],
            "source": "workflow_scope_single",
            "error": "",
        }

    return {"mode": None, "fwd_files": [], "rev_files": [], "single_files": [], "source": None, "error": "No valid DADA2 inputs found."}


def _resolve_boldigger_python():
    candidate = BOLDIGGER_PYTHON or sys.executable
    if candidate and (os.path.exists(candidate) or shutil.which(candidate)):
        return candidate
    return sys.executable


def _check_python_module(python_bin, module_name):
    try:
        result = subprocess.run(
            [python_bin, "-c", f"import {module_name}"],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except Exception as exc:
        return False, str(exc)
    if result.returncode != 0:
        return False, (result.stderr or result.stdout or f"Unable to import {module_name}").strip()
    return True, ""


def _blast_version_info():
    version = _safe_version([BLASTN_PATH, "-version"])
    return version or "Unknown BLAST+ version"


def _run_blastn_with_retry(query_fasta, db_target, output_path, max_hits=5, evalue="1e-10", remote=False, timeout_sec=BLAST_TIMEOUT_SEC, retries=3, backoff_sec=2):
    attempts = []
    last_stdout = ""
    last_stderr = ""
    last_code = 1
    total_attempts = retries if remote else 1
    for attempt in range(1, total_attempts + 1):
        stdout, stderr, code = run_blastn(query_fasta, db_target, output_path, max_hits=max_hits, evalue=evalue, remote=remote, timeout_sec=timeout_sec)
        attempts.append({"attempt": attempt, "returncode": code, "stderr": stderr, "stdout": stdout})
        last_stdout, last_stderr, last_code = stdout, stderr, code
        if code == 0 and os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            return stdout, stderr, code, attempts
        if not remote or attempt >= total_attempts:
            break
        time.sleep(backoff_sec * (2 ** (attempt - 1)))
    return last_stdout, last_stderr, last_code, attempts


def _resolve_output_path(project_state, step_name, output_type, filename=None, sample_id=None):
    base = _get_active_root(project_state)
    if base == UPLOAD_DIRECTORY:
        base = os.path.abspath(base)

    steps_root = os.path.join(base, "steps")
    results_root = os.path.join(base, "results")
    os.makedirs(steps_root, exist_ok=True)
    os.makedirs(results_root, exist_ok=True)

    if output_type == "intermediate":
        target_dir = os.path.join(steps_root, step_name)
        if sample_id:
            target_dir = os.path.join(target_dir, _safe_name(sample_id))
        os.makedirs(target_dir, exist_ok=True)
        return os.path.join(target_dir, filename) if filename else target_dir

    if output_type == "final":
        target_dir = os.path.join(steps_root, step_name, "outputs") if step_name else results_root
        if sample_id:
            target_dir = os.path.join(target_dir, _safe_name(sample_id))
        os.makedirs(target_dir, exist_ok=True)
        return os.path.join(target_dir, filename) if filename else target_dir

    raise ValueError(f"Unknown output_type: {output_type}")


def _resolve_step_input(candidate_paths):
    for source_label, path in candidate_paths or []:
        if path and os.path.exists(path) and os.path.getsize(path) > 0:
            return path, source_label, {
                "source": source_label,
                "path": path,
                "resolved_at": datetime.now().isoformat(timespec="seconds"),
            }
    return None, None, None


def _resolve_preprocess_sources(uploaded_files, samples_store, active_sample):
    data = uploaded_files or {}
    sample_input = _resolve_sample_input(samples_store, active_sample)
    se_uploads = _normalize_uploaded_paths(data.get("preprocess", []))
    se_input = se_uploads[0] if se_uploads else sample_input.get("se")
    r1_input = data.get("preprocess_r1") or sample_input.get("r1")
    r2_input = data.get("preprocess_r2") or sample_input.get("r2")
    return {
        "se": se_input,
        "r1": r1_input,
        "r2": r2_input,
    }


def _resolve_sample_analysis_source(sample_id, sample_info, output_dir, project_state=None):
    info = sample_info or {}
    outputs = info.get("outputs", {}) if isinstance(info, dict) else {}
    paths = info.get("paths", {}) if isinstance(info, dict) else {}

    preprocess_base = _step_base(project_state, "preprocess") if project_state else ""
    pear_merged = ""

    # First priority: sample-specific PEAR merged file
    if preprocess_base and isinstance(sample_id, str):
        safe_id = _safe_name(sample_id)
        import glob as _glob
        pear_pattern = os.path.join(preprocess_base, "**", "merged.assembled.fastq")

        for pear_candidate in _glob.glob(pear_pattern, recursive=True):
            parent = os.path.basename(os.path.dirname(pear_candidate)).lower()
            members = info.get("members") if isinstance(info, dict) else []
            if safe_id.lower() in parent or any(
                isinstance(m, str) and m.lower().replace(" ", "_") in parent
                for m in (members or [])
            ):
                pear_merged = pear_candidate
                break

        if not pear_merged:
            all_merged = _glob.glob(pear_pattern, recursive=True)
            if all_merged:
                pear_merged = all_merged[0]

    # Fallbacks after PEAR
    candidate_paths = [
        pear_merged if pear_merged else None,
        outputs.get("input_fasta"),
        outputs.get("asv_fasta"),
        paths.get("se_fasta"),
        paths.get("se"),
        paths.get("r1"),
        paths.get("r2"),
    ]

    source_path = _resolve_local_source(candidate_paths)
    if not source_path:
        return None

    lower = source_path.lower()
    if lower.endswith((".fasta", ".fa")):
        return source_path

    if lower.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
        temp_dir = os.path.join(output_dir, "pooled_inputs")
        os.makedirs(temp_dir, exist_ok=True)
        temp_fasta = os.path.join(temp_dir, f"{_safe_name(sample_id)}.fasta")
        if fastq_to_fasta(source_path, temp_fasta):
            return temp_fasta

    return None


def _build_pooled_fasta_for_analysis(samples_store, workflow_scope, output_dir, project_state=None):
    os.makedirs(output_dir, exist_ok=True)
    combined_fasta = os.path.join(output_dir, "pooled_sequences.fasta")
    sample_counts = {}

    from Bio import SeqIO

    with open(combined_fasta, "w", encoding="utf-8") as out_handle:
        for sample_id, sample_info in _iter_workflow_samples(samples_store, workflow_scope):
            source_path = _resolve_sample_analysis_source(sample_id, sample_info, output_dir, project_state)
            if not source_path or not os.path.exists(source_path) or os.path.getsize(source_path) == 0:
                continue

            seq_count = 0
            for record in SeqIO.parse(source_path, "fasta"):
                record.id = f"{_safe_name(sample_id)}|{record.id}"
                record.name = record.id
                record.description = ""
                SeqIO.write(record, out_handle, "fasta")
                seq_count += 1
            if seq_count > 0:
                sample_counts[sample_id] = seq_count

    if not os.path.exists(combined_fasta) or os.path.getsize(combined_fasta) == 0:
        return None, None

    return combined_fasta, sample_counts


def _normalize_samples_store(samples_store):
    """Merge/normalize an existing samples_store so PE pairs become single logical samples.

    This is a migration helper to coalesce entries that represent R1/R2 members
    into one sample keyed by a pair_key.
    """
    if not samples_store or not isinstance(samples_store, dict):
        return samples_store or {}

    merged = {}
    seen_pair = {}

    def _pair_key_from_info(info):
        name = info.get('display_name') or info.get('stem') or ''
        return _infer_marker_from_name(name) if False else re.sub(r"[^A-Za-z0-9]+", "_", re.sub(r"(?i)(?:[_\.-]?(?:raw[_\.-]?)?(?:r|fwd|rev)?[ _\.-]?(?:1|2))$", "", os.path.splitext(str(name))[0])).strip("_")

    for sid, info in list(samples_store.items()):
        pk = info.get('pair_key') or _pair_key_from_info(info)
        if not pk:
            pk = sid
        if pk in seen_pair:
            target = seen_pair[pk]
            # merge members
            merged[target].setdefault('members', [])
            for m in info.get('members', []) if isinstance(info.get('members', []), list) else [info.get('display_name') or info.get('stem') or sid]:
                if m and m not in merged[target]['members']:
                    merged[target]['members'].append(m)
            # merge paths
            for pkey, pval in (info.get('paths') or {}).items():
                if pval:
                    merged[target]['paths'][pkey] = pval
            # prefer existing kind PE if any
            if merged[target].get('kind') != 'PE' and info.get('kind') == 'PE':
                merged[target]['kind'] = 'PE'
                merged[target]['mode'] = 'PE'
            if not merged[target].get('internal_sample_id'):
                merged[target]['internal_sample_id'] = info.get('internal_sample_id') or pk
            if not merged[target].get('display_label'):
                merged[target]['display_label'] = merged[target].get('display_name') or merged[target].get('pair_key') or target
        else:
            # create new merged entry
            new_sid = sid
            merged[new_sid] = {
                'display_name': info.get('display_name') or info.get('stem') or new_sid,
                'display_label': info.get('display_label') or info.get('display_name') or info.get('stem') or new_sid,
                'stem': info.get('stem') or new_sid,
                'pair_key': pk,
                'internal_sample_id': info.get('internal_sample_id') or pk,
                'kind': info.get('kind', 'SE'),
                'mode': info.get('mode', info.get('kind', 'SE')),
                'marker': info.get('marker'),
                'members': list(info.get('members', [info.get('display_name') or info.get('stem') or new_sid])) if info.get('members') else [info.get('display_name') or info.get('stem') or new_sid],
                'paths': dict(info.get('paths', {}) or {}),
                'outputs': dict(info.get('outputs', {}) or {}),
                'meta': dict(info.get('meta', {}) or {})
            }
            seen_pair[pk] = new_sid

    return merged


def _parse_vsearch_uc_to_matrix(uc_file, sample_counts):
    if not uc_file or not os.path.exists(uc_file) or not sample_counts:
        return None

    otu_counts = {}
    otu_order = []

    with open(uc_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            otu_key = parts[1]
            query_label = parts[8]
            if "|" not in query_label:
                continue
            sample_id = query_label.split("|", 1)[0]
            if sample_id not in sample_counts:
                continue
            if otu_key not in otu_counts:
                otu_counts[otu_key] = {}
                otu_order.append(otu_key)
            otu_counts[otu_key][sample_id] = otu_counts[otu_key].get(sample_id, 0) + 1

    if not otu_counts:
        return None

    rows = []
    for otu_key in otu_order:
        row = {sample_id: otu_counts[otu_key].get(sample_id, 0) for sample_id in sample_counts.keys()}
        row["OTU"] = f"OTU_{len(rows) + 1}"
        rows.append(row)

    matrix = pd.DataFrame(rows)
    if matrix.empty:
        return None
    matrix = matrix.set_index("OTU")
    return matrix


def _normalize_uploaded_paths(filenames):
    if not filenames:
        return []
    if isinstance(filenames, str):
        filenames = [filenames]
    if not isinstance(filenames, list):
        return []
    normalized = []
    for item in filenames:
        if item:
            normalized.append(str(item))
    return normalized


def _as_selection_list(value):
    if value in (None, "", []):
        return []
    if isinstance(value, list):
        return [str(v) for v in value if v]
    return [str(value)]


def _normalize_active_sample_value(value):
    selected = _as_selection_list(value)
    if not selected or "__ALL__" in selected:
        return "__ALL__"
    if len(selected) == 1:
        return selected[0]
    return selected


def _ctx_uploaded_files_store():
    try:
        current = (dash.callback_context.states or {}).get("uploaded-files.data") or {}
    except Exception:
        current = {}
    if not isinstance(current, dict):
        return {}
    return current.copy()


def _merge_uploaded_files_from_ctx(patch):
    merged = _ctx_uploaded_files_store()
    if not isinstance(patch, dict):
        return merged
    merged.update(patch)
    return merged


def _project_session_path(project_dir):
    return os.path.join(project_dir, "session_state.json")


def _run_session_path(run_dir):
    return os.path.join(run_dir, "session.json")


STEP_STATUS_PENDING = "pending"
STEP_STATUS_RUNNING = "running"
STEP_STATUS_SUCCESS = "success"
STEP_STATUS_FAILED = "failed"
STEP_STATUS_SKIPPED = "skipped"


def _step_manifest_path(project_state, step_name):
    return _resolve_output_path(project_state, step_name, "final", f"{step_name}_manifest.json")


def _final_manifest_path(project_state):
    return _resolve_output_path(project_state, "analysis", "final", "final_manifest.json")


def _path_is_valid(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _step_output_entry(path, label=None):
    return {
        "label": label or os.path.basename(path) if path else "",
        "path": path or "",
        "exists": bool(path) and os.path.exists(path),
        "size": int(os.path.getsize(path)) if path and os.path.exists(path) else 0,
        "valid": _path_is_valid(path),
    }


def _build_step_manifest(step_name, status, outputs=None, inputs=None, error=None, notes=None):
    return {
        "step": step_name,
        "status": status,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "outputs": outputs or {},
        "inputs": inputs or {},
        "error": error or "",
        "notes": notes or "",
    }


def _write_step_manifest(project_state, step_name, manifest):
    if not isinstance(manifest, dict):
        return False
    return _write_json(_step_manifest_path(project_state, step_name), manifest)


def _read_step_manifest(project_state, step_name, default=None):
    return _read_json(_step_manifest_path(project_state, step_name), default=default)


def _manifest_path_value(payload, *keys):
    current = payload or {}
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current if _path_is_valid(current) else None


def _build_run_session_metadata(project_state, run_config, workflow_scope, samples_store, completed_at=None, outputs=None):
    state = project_state or {}
    workflow_samples = _get_workflow_scope_samples(workflow_scope)
    if not workflow_samples:
        workflow_samples = list((samples_store or {}).keys())

    metadata = {
        "run_id": os.path.basename(state.get("run_dir") or ""),
        "project_name": state.get("project_name", ""),
        "project_dir": state.get("project_dir", ""),
        "run_dir": state.get("run_dir", ""),
        "created_at": (run_config or {}).get("created_at") or datetime.now().isoformat(timespec="seconds"),
        "completed_at": completed_at or "",
        "samples_selected": workflow_samples,
        "steps_executed": (run_config or {}).get("steps", []),
        "parameters": run_config or {},
        "tool_versions": _collect_tool_versions(),
        "outputs": outputs or {},
        "execution_log": [],
    }
    return metadata


def _write_run_session_metadata(project_state, run_config, workflow_scope, samples_store, completed_at=None, outputs=None):
    state = project_state or {}
    run_dir = state.get("run_dir")
    if not run_dir:
        return False
    metadata = _build_run_session_metadata(state, run_config, workflow_scope, samples_store, completed_at=completed_at, outputs=outputs)
    try:
        _write_json(_run_session_path(run_dir), metadata)
        return True
    except Exception:
        return False


def _update_run_session_metadata(project_state, updates):
    state = project_state or {}
    run_dir = state.get("run_dir")
    if not run_dir:
        return False

    session_json_path = _run_session_path(run_dir)
    session_payload = _read_json(session_json_path, default={}) or {}
    if not isinstance(session_payload, dict):
        session_payload = {}

    if not isinstance(updates, dict) or not updates:
        return False

    for key, value in updates.items():
        session_payload[key] = value

    try:
        _write_json(session_json_path, session_payload)
        return True
    except Exception:
        return False


def _list_existing_projects():
    projects = []
    if not os.path.exists(PROJECTS_ROOT):
        return projects
    for entry in sorted(os.listdir(PROJECTS_ROOT)):
        pdir = os.path.join(PROJECTS_ROOT, entry)
        if not os.path.isdir(pdir):
            continue
        meta_path = os.path.join(pdir, "project.json")
        if os.path.exists(meta_path):
            projects.append({"label": entry, "value": entry})
    return projects


def _latest_run_dir(project_dir):
    runs_dir = os.path.join(project_dir, "runs")
    if not os.path.exists(runs_dir):
        return None
    candidates = []
    for entry in os.listdir(runs_dir):
        rdir = os.path.join(runs_dir, entry)
        if os.path.isdir(rdir):
            candidates.append((os.path.getmtime(rdir), rdir))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def _read_json(path, default=None):
    if not os.path.exists(path):
        return default
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default


def _write_json(path, payload):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(prefix="tmp_session_", suffix=".json", dir=os.path.dirname(path))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def _count_fasta_seqs(path):
    if not path or not os.path.exists(path):
        return 0
    count = 0
    try:
        with open(path, 'r', encoding='utf-8') as fh:
            for line in fh:
                if line.startswith('>'):
                    count += 1
    except Exception:
        return 0
    return count


def _get_directory_size(path):
    if not path or not os.path.exists(path):
        return 0
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            file_path = os.path.join(root, name)
            try:
                total += os.path.getsize(file_path)
            except OSError:
                continue
    return total


def _apply_cleanup_policy(run_dir, delete_after_steps=None):
    if KEEP_INTERMEDIATES or not run_dir:
        return []
    steps_to_delete = [step for step in (delete_after_steps if delete_after_steps is not None else DELETE_AFTER_STEPS) if step]
    deleted = []
    for step_name in steps_to_delete:
        step_dir = os.path.join(run_dir, "steps", step_name)
        if os.path.exists(step_dir):
            shutil.rmtree(step_dir, ignore_errors=True)
            deleted.append(step_name)
    return deleted


def _save_project_session(project_state, run_config, uploaded_files, active_sample, samples_store, workflow_scope):
    state = project_state or {}
    project_dir = state.get("project_dir")
    if not project_dir:
        return False
    payload = {
        "saved_at": datetime.now().isoformat(timespec="seconds"),
        "project_state": state,
        "run_config": run_config or {},
        "uploaded_files": uploaded_files or {"preprocess_mode": "se"},
        "active_sample": active_sample if active_sample is not None else "__ALL__",
        "samples_store": samples_store or {},
        "workflow_scope": workflow_scope or {}
    }
    _write_json(_project_session_path(project_dir), payload)
    return True


def _load_project_session(project_dir):
    payload = _read_json(_project_session_path(project_dir), default={}) or {}
    state = payload.get("project_state", {}) if isinstance(payload, dict) else {}
    state = state if isinstance(state, dict) else {}
    run_config = payload.get("run_config", {}) if isinstance(payload, dict) else {}
    run_config = run_config if isinstance(run_config, dict) else {}
    uploaded_files = payload.get("uploaded_files", {"preprocess_mode": "se"}) if isinstance(payload, dict) else {"preprocess_mode": "se"}
    uploaded_files = uploaded_files if isinstance(uploaded_files, dict) else {"preprocess_mode": "se"}
    active_sample = payload.get("active_sample", "__ALL__") if isinstance(payload, dict) else "__ALL__"
    samples_store = payload.get("samples_store", {}) if isinstance(payload, dict) else {}
    samples_store = samples_store if isinstance(samples_store, dict) else {}
    workflow_scope = payload.get("workflow_scope", {}) if isinstance(payload, dict) else {}
    workflow_scope = workflow_scope if isinstance(workflow_scope, dict) else {}
    return state, run_config, uploaded_files, active_sample, samples_store, workflow_scope


def _iter_target_samples(samples_store, active_sample):
    samples = samples_store or {}
    if not isinstance(samples, dict) or not samples:
        return []

    if isinstance(active_sample, list):
        selected = [sid for sid in active_sample if sid in samples]
        return [(sid, samples.get(sid) or {}) for sid in selected]

    if active_sample and active_sample != "__ALL__" and active_sample in samples:
        return [(active_sample, samples.get(active_sample) or {})]
    return list(samples.items())


def _get_workflow_scope_samples(workflow_scope):
    """Extract sample IDs from workflow scope store.
    
    Args:
        workflow_scope: dict with structure {"samples_selected": [...], "selected_at": "...", "pe_pairs": {...}}
    
    Returns:
        list of sample IDs, or empty list if "__ALL__" or None
    """
    scope = workflow_scope or {}
    samples = scope.get("samples_selected", [])
    if not samples or "__ALL__" in samples:
        return []  # Will be handled by caller (iterate all)
    return samples


def _group_samples_for_options(samples_store):
    """Return a flattened list of options for selectors that groups samples by
    biological id + marker + mode but preserves distinct technical pairs as
    separate selectable values. Each option value remains the original sample SID.
    """
    samples = samples_store or {}
    groups = {}
    options = []

    for sid, info in samples.items():
        display = info.get('display_label') or info.get('display_name') or info.get('pair_key') or info.get('stem') or sid
        marker = info.get('marker') or _infer_marker_from_name(display)
        kind = info.get('kind', info.get('mode', 'SE'))

        bio_key = f"{display}|{kind}|{sid}"
        groups.setdefault(bio_key, []).append((sid, info))

    for bio_key, members in sorted(groups.items(), key=lambda x: x[0]):
        # members are distinct technical sets (e.g. multiple runs/pairs)
        display_name, kind, _ = bio_key.split("|")
        prefix = display_name
        for idx, (sid, info) in enumerate(sorted(members, key=lambda x: x[0])):
            label = prefix
            if kind == 'PE' or info.get('mode') == 'PE':
                if len(members) > 1:
                    label = f"{label} [set {idx + 1}]"
            options.append({"label": label, "value": sid})
    return options


def _infer_marker_from_name(name):
    if not name:
        return None
    n = str(name).lower()
    if any(x in n for x in ["coi", "cox1", "mtcoi"]):
        return "COI"
    if any(x in n for x in ["16s", "rrna", "ssu", "16s_r", "16s_rna"]):
        return "16S"
    if any(x in n for x in ["18s", "18s_r", "euk", "ssu"]):
        return "18S"
    if "its" in n:
        return "ITS"
    return None


def _infer_marker_from_samples(samples_store, workflow_scope):
    scope = workflow_scope or {}
    selected = scope.get('samples_selected', [])
    samples = samples_store or {}
    candidates = []
    if not selected or "__ALL__" in selected:
        # consider all
        for sid, info in samples.items():
            candidates.append(info.get('display_name') or info.get('stem') or sid)
    else:
        for sid in selected:
            info = samples.get(sid, {})
            candidates.append(info.get('display_name') or info.get('stem') or sid)

    # Infer by majority of per-sample inferences
    inferred = []
    for c in candidates:
        m = _infer_marker_from_name(c)
        if m:
            inferred.append(m)
    if not inferred:
        return None
    return max(set(inferred), key=inferred.count)


def _iter_workflow_samples(samples_store, workflow_scope):
    """Iterate target samples using workflow scope (global run-wide selection).
    
    Args:
        samples_store: Complete dict of all samples
        workflow_scope: dict with structure {"samples_selected": [...], ...}
    
    Returns:
        list of (sample_id, sample_info) tuples for samples in workflow scope
    """
    scope = workflow_scope or {}
    selected_ids = scope.get("samples_selected", [])
    
    if not selected_ids or "__ALL__" in selected_ids:
        # Return all samples from store
        return list((samples_store or {}).items())
    
    # Return only selected samples
    result = []
    for sid in selected_ids:
        if sid in (samples_store or {}):
            result.append((sid, samples_store[sid]))
    return result


def _create_session_results_zip(project_state):
    state = project_state or {}
    run_dir = state.get("run_dir")
    project_dir = state.get("project_dir")

    if not run_dir or not os.path.isdir(run_dir):
        _set_export_error("Export not created: no completed run directory was available.")
        return None

    exports_dir = _project_exports_dir(project_dir)
    if not exports_dir:
        _set_export_error("Export not created: could not prepare the project exports directory.")
        return None

    zip_name = f"PARAFISH_{_safe_name(state.get('project_name', 'Project'))}_{_safe_name(os.path.basename(run_dir))}.zip"
    zip_path = os.path.join(exports_dir, zip_name)

    try:
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
            for root_name in ("steps", "results"):
                root_dir = os.path.join(run_dir, root_name)
                if not os.path.isdir(root_dir):
                    continue
                for current_root, dirnames, filenames in os.walk(root_dir):
                    dirnames[:] = [name for name in dirnames if not _is_export_blocked_component(name)]
                    for filename in filenames:
                        file_path = os.path.join(current_root, filename)
                        if not _should_include_export_file(file_path, run_dir):
                            continue
                        archive.write(file_path, os.path.relpath(file_path, run_dir))
    except OSError as exc:
        try:
            if os.path.exists(zip_path):
                os.remove(zip_path)
        except OSError:
            pass
        _set_export_error(f"Export not created: unable to write ZIP archive ({exc}).")
        return None

    if not _path_is_valid(zip_path):
        try:
            if os.path.exists(zip_path):
                os.remove(zip_path)
        except OSError:
            pass
        _set_export_error("Export not created: no eligible step/results files were found.")
        return None

    _set_export_error("")
    StorageManager.cleanup_export_archives(project_dir, run_dir)
    return zip_path

def _resolve_sample_input(samples_store, active_sample, prefer_pe=True):
    if active_sample in (None, "__ALL__"):
        return {}
    if isinstance(active_sample, list):
        active_sample = active_sample[0] if active_sample else None
        if not active_sample:
            return {}
    if not samples_store or not active_sample:
        return {}
    sample = samples_store.get(active_sample, {})
    paths = sample.get("paths", {})
    if prefer_pe and paths.get("r1") and paths.get("r2"):
        return {"r1": paths.get("r1"), "r2": paths.get("r2")}
    if paths.get("se"):
        return {"se": paths.get("se")}
    return {}


def _sample_fastq_candidates(sample_info):
    info = sample_info or {}
    paths = info.get("paths", {}) if isinstance(info, dict) else {}

    # Prefer explicit single-end input when present; otherwise use paired reads.
    if paths.get("se"):
        return [paths.get("se")]

    candidates = []
    if paths.get("r1"):
        candidates.append(paths.get("r1"))
    if paths.get("r2") and paths.get("r2") != paths.get("r1"):
        candidates.append(paths.get("r2"))
    return candidates


def _build_query_fasta_from_samples(samples_store, active_sample, output_dir):
    targets = _iter_target_samples(samples_store, active_sample)
    if not targets:
        return None, "No samples are available to build a query FASTA."

    os.makedirs(output_dir, exist_ok=True)

    fastq_files = []
    seen = set()
    for _, sample_info in targets:
        for fp in _sample_fastq_candidates(sample_info):
            if fp and fp not in seen:
                seen.add(fp)
                fastq_files.append(fp)

    if not fastq_files:
        return None, "No FASTQ files were found in selected project samples for DB classification."

    converted_fastas = []
    for idx, fastq_path in enumerate(fastq_files, start=1):
        if not os.path.exists(fastq_path) or os.path.getsize(fastq_path) == 0:
            continue
        stem = os.path.basename(fastq_path)
        stem = stem.replace(".fastq.gz", "").replace(".fq.gz", "").replace(".fastq", "").replace(".fq", "")
        out_fasta = os.path.join(output_dir, f"sample_{idx:03d}_{_safe_name(stem)}.fasta")
        if fastq_to_fasta(fastq_path, out_fasta):
            converted_fastas.append(out_fasta)

    if not converted_fastas:
        return None, "Could not convert project FASTQ files to FASTA for DB classification."

    if len(converted_fastas) == 1:
        return converted_fastas[0], None

    combined_fasta = os.path.join(output_dir, "query_from_samples_combined.fasta")
    if not write_combined_fasta(converted_fastas, combined_fasta):
        return None, "Failed to combine converted FASTA files for DB classification."
    return combined_fasta, None


def _upload_callback_payload(filenames, store_key, label, progress_done_key):
    normalized_files = _normalize_uploaded_paths(filenames)
    uploaded_files_update = _merge_uploaded_files_from_ctx({})

    if not normalized_files:
        uploaded_files_update[store_key] = []
        return html.Div(f"No files uploaded for {label}."), dash.no_update, uploaded_files_update

    uploaded_files_update[store_key] = [normalized_files[0]]
    progress_update = {
        'qc_done': True,
        'preprocess_done': True,
        'analysis_done': True,
        'db_done': True,
        'phylo_done': True,
        'results_done': False
    }
    if progress_done_key in progress_update:
        progress_update[progress_done_key] = True

    status = html.Ul([html.Li(f"{label}: {os.path.basename(normalized_files[0])}")])
    return status, progress_update, uploaded_files_update

def read_file_preview(path, n_lines=30):
    if not os.path.exists(path):
        return "File not found."
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            lines = []
            for i, line in enumerate(f):
                if i >= n_lines:
                    lines.append("... (truncated output) ...")
                    break
                lines.append(line.rstrip())
        return "\n".join(lines)
    except Exception as e:
        return f"Error reading file: {str(e)}"

def parse_fastqc_basic_stats(fastqc_data_path):
    stats = {}
    in_basic = False
    try:
        with open(fastqc_data_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">>Basic Statistics"):
                    in_basic = True
                    continue
                if in_basic and line.startswith(">>"):
                    break
                if in_basic and line and not line.startswith("#"):
                    parts = line.split("\t")
                    if len(parts) == 2:
                        stats[parts[0]] = parts[1]
    except Exception:
        return {}
    return stats


def _fastqc_report_paths(input_file, output_dir):
    basename = os.path.basename(input_file)
    if basename.endswith(".fastq.gz"):
        stem = basename[:-9]
    elif basename.endswith(".fastq"):
        stem = basename[:-6]
    elif basename.endswith(".fq.gz"):
        stem = basename[:-6]
    elif basename.endswith(".fq"):
        stem = basename[:-3]
    else:
        stem = os.path.splitext(basename)[0]
    return {
        "html": os.path.join(output_dir, f"{stem}_fastqc.html"),
        "zip": os.path.join(output_dir, f"{stem}_fastqc.zip"),
    }


def _parse_fastqc_basic_stats_lines(lines):
    stats = {}
    in_basic = False
    for raw_line in lines:
        line = raw_line.decode("utf-8", errors="replace") if isinstance(raw_line, (bytes, bytearray)) else str(raw_line)
        line = line.strip()
        if line.startswith(">>Basic Statistics"):
            in_basic = True
            continue
        if in_basic and line.startswith(">>"):
            break
        if in_basic and line and not line.startswith("#"):
            parts = line.split("\t")
            if len(parts) == 2:
                stats[parts[0]] = parts[1]
    return stats


def _parse_fastqc_basic_stats_from_zip(zip_path):
    if not zip_path or not os.path.exists(zip_path) or os.path.getsize(zip_path) == 0:
        return {}
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            candidates = [name for name in zf.namelist() if name.endswith("fastqc_data.txt")]
            if not candidates:
                return {}
            with zf.open(candidates[0], "r") as handle:
                return _parse_fastqc_basic_stats_lines(handle)
    except Exception:
        return {}

def _download_file(url, dest_path):
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    urllib.request.urlretrieve(url, dest_path)
    return dest_path

def _gunzip_file(gz_path, out_path):
    with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


def _file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _bold_cache_key(query_fasta_path, db_nr, mode, thresholds):
    thresholds_token = ""
    if thresholds and str(thresholds).strip():
        thresholds_token = " ".join(str(thresholds).strip().split())
    payload = "|".join([
        _file_sha256(query_fasta_path),
        str(db_nr),
        str(mode),
        thresholds_token,
    ])
    return hashlib.sha256(payload.encode("utf-8", errors="replace")).hexdigest()


def _resolve_bold_cache_hit(cache_key):
    cache_dir = os.path.join(BOLD_CACHE_ROOT, cache_key)
    xlsx_path = os.path.join(cache_dir, "result.xlsx")
    parquet_path = os.path.join(cache_dir, "result.parquet")
    meta_path = os.path.join(cache_dir, "meta.json")
    has_xlsx = os.path.exists(xlsx_path) and os.path.getsize(xlsx_path) > 0
    has_parquet = os.path.exists(parquet_path) and os.path.getsize(parquet_path) > 0
    if not has_xlsx and not has_parquet:
        return None
    return {
        "cache_dir": cache_dir,
        "xlsx": xlsx_path if has_xlsx else None,
        "parquet": parquet_path if has_parquet else None,
        "meta": meta_path if os.path.exists(meta_path) else None,
    }


def _save_bold_cache(cache_key, source_xlsx=None, source_parquet=None, metadata=None):
    cache_dir = os.path.join(BOLD_CACHE_ROOT, cache_key)
    os.makedirs(cache_dir, exist_ok=True)
    if source_xlsx and os.path.exists(source_xlsx) and os.path.getsize(source_xlsx) > 0:
        shutil.copy2(source_xlsx, os.path.join(cache_dir, "result.xlsx"))
    if source_parquet and os.path.exists(source_parquet) and os.path.getsize(source_parquet) > 0:
        shutil.copy2(source_parquet, os.path.join(cache_dir, "result.parquet"))
    if metadata:
        with open(os.path.join(cache_dir, "meta.json"), "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)
    return cache_dir

def _candidate_db_source_paths(db_name, dest_dir, explicit_local_path=None):
    suffixes = (".fasta.gz", ".fa.gz", ".fasta", ".fa")
    candidates = []

    for path in [
        explicit_local_path,
        os.path.join(dest_dir, f"{db_name}.fasta.gz"),
        os.path.join(dest_dir, f"{db_name}.fa.gz"),
        os.path.join(dest_dir, f"{db_name}.fasta"),
        os.path.join(dest_dir, f"{db_name}.fa"),
    ]:
        if path:
            candidates.append(path)

    for root_dir in [DB_ROOT]:
        if not os.path.isdir(root_dir):
            continue
        for entry in sorted(os.listdir(root_dir)):
            entry_lower = entry.lower()
            if db_name in entry_lower and entry_lower.endswith(suffixes):
                candidates.append(os.path.join(root_dir, entry))

    unique_candidates = []
    seen = set()
    for path in candidates:
        normalized = os.path.abspath(path)
        if normalized in seen:
            continue
        seen.add(normalized)
        unique_candidates.append(path)
    return unique_candidates


def ensure_fasta_db(db_name, url, dest_dir, local_path=None):
    fasta_path = os.path.join(dest_dir, f"{db_name}.fasta")
    if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 0:
        return fasta_path, None
    os.makedirs(dest_dir, exist_ok=True)
    gz_path = os.path.join(dest_dir, f"{db_name}.fasta.gz")

    for source_path in _candidate_db_source_paths(db_name, dest_dir, explicit_local_path=local_path):
        if not os.path.exists(source_path) or os.path.getsize(source_path) == 0:
            continue
        try:
            if source_path.lower().endswith(".gz"):
                source_gz = source_path
                if os.path.abspath(source_path) != os.path.abspath(gz_path):
                    shutil.copy2(source_path, gz_path)
                    source_gz = gz_path
                _gunzip_file(source_gz, fasta_path)
            else:
                if os.path.abspath(source_path) != os.path.abspath(fasta_path):
                    shutil.copy2(source_path, fasta_path)
            if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 0:
                return fasta_path, None
        except Exception as e:
            return None, f"Failed to prepare local {db_name} database from {source_path}: {str(e)}"

    if not url:
        return None, f"Missing local file or download URL for {db_name}."

    if not ALLOW_DB_AUTO_DOWNLOAD:
        return None, (
            f"{db_name.upper()} is not available locally under global DB root ({DB_ROOT}) and auto-download is disabled. "
            f"Provide local file via env var/path or enable ALLOW_DB_AUTO_DOWNLOAD=1."
        )

    try:
        _download_file(url, gz_path)
        _gunzip_file(gz_path, fasta_path)
        return fasta_path, None
    except Exception as e:
        return None, f"Failed to download {db_name}: {str(e)}"


def _resolve_local_source(candidates):
    for path in candidates:
        if path and os.path.exists(path) and os.path.getsize(path) > 0:
            return path
    return None


def _prepare_local_fasta_or_gz(source_path, dest_fasta, dest_gz=None):
    if source_path.lower().endswith(".gz"):
        src_gz = source_path
        if dest_gz:
            os.makedirs(os.path.dirname(dest_gz), exist_ok=True)
            if os.path.abspath(source_path) != os.path.abspath(dest_gz):
                shutil.copy2(source_path, dest_gz)
                src_gz = dest_gz
        _gunzip_file(src_gz, dest_fasta)
    else:
        if os.path.abspath(source_path) != os.path.abspath(dest_fasta):
            shutil.copy2(source_path, dest_fasta)


def _normalize_pr2_mothur_db(fasta_in, tax_in, out_fasta, out_map):
    tax_by_orig = {}
    with open(tax_in, "r", encoding="utf-8", errors="replace") as tf:
        for line in tf:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                tax_by_orig[parts[0]] = parts[1]

    def _fasta_records(path):
        header = None
        seq_parts = []
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_parts)
                    header = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if header is not None:
                yield header, "".join(seq_parts)

    count = 0
    with open(out_fasta, "w", encoding="utf-8") as fo, open(out_map, "w", encoding="utf-8") as fm:
        fm.write("normalized_id\toriginal_id\ttaxonomy\n")
        for orig_id, seq in _fasta_records(fasta_in):
            if not seq:
                continue
            count += 1
            norm_id = f"PR2_{count:06d}"
            fo.write(f">{norm_id}\n{seq}\n")
            fm.write(f"{norm_id}\t{orig_id}\t{tax_by_orig.get(orig_id, 'NA')}\n")
    return count


def prepare_pr2_mothur_db(pr2_dir):
    os.makedirs(pr2_dir, exist_ok=True)

    fasta_candidates = [
        DEFAULT_PR2_MOTHUR_FASTA_PATH,
        os.path.join(pr2_dir, "pr2_version_5.1.1_SSU_mothur.fasta.gz"),
        os.path.join(pr2_dir, "pr2_version_5.1.1_SSU_mothur.fasta"),
        os.path.join(os.getcwd(), "pr2_version_5.1.1_SSU_mothur.fasta.gz"),
        os.path.join(os.getcwd(), "pr2_version_5.1.1_SSU_mothur.fasta"),
    ]
    tax_candidates = [
        DEFAULT_PR2_MOTHUR_TAX_PATH,
        os.path.join(pr2_dir, "pr2_version_5.1.1_SSU_mothur.tax.gz"),
        os.path.join(pr2_dir, "pr2_version_5.1.1_SSU_mothur.tax"),
        os.path.join(os.getcwd(), "pr2_version_5.1.1_SSU_mothur.tax.gz"),
        os.path.join(os.getcwd(), "pr2_version_5.1.1_SSU_mothur.tax"),
    ]

    mothur_src = _resolve_local_source(fasta_candidates)
    tax_src = _resolve_local_source(tax_candidates)

    if not mothur_src:
        return None, None, "PR2 mothur FASTA not found. Place pr2_version_5.1.1_SSU_mothur.fasta.gz in project root or databases/pr2/5_1_1/."
    if not tax_src:
        return None, None, "PR2 mothur TAX not found. Place pr2_version_5.1.1_SSU_mothur.tax.gz in project root or databases/pr2/5_1_1/."

    mothur_fasta = os.path.join(pr2_dir, "pr2_SSU_mothur.fasta")
    mothur_tax = os.path.join(pr2_dir, "pr2_SSU_mothur.tax")

    try:
        _prepare_local_fasta_or_gz(mothur_src, mothur_fasta, dest_gz=os.path.join(pr2_dir, "pr2_SSU_mothur.fasta.gz"))
        _prepare_local_fasta_or_gz(tax_src, mothur_tax, dest_gz=os.path.join(pr2_dir, "pr2_SSU_mothur.tax.gz"))
    except Exception as e:
        return None, None, f"Failed to prepare PR2 mothur files: {e}"

    normalized_fasta = os.path.join(pr2_dir, "pr2.normal.fa")
    taxmap_tsv = os.path.join(pr2_dir, "pr2_taxmap.tsv")

    if os.path.exists(normalized_fasta) and os.path.exists(taxmap_tsv) and os.path.getsize(normalized_fasta) > 0 and os.path.getsize(taxmap_tsv) > 0:
        return normalized_fasta, taxmap_tsv, None

    try:
        n_records = _normalize_pr2_mothur_db(mothur_fasta, mothur_tax, normalized_fasta, taxmap_tsv)
    except Exception as e:
        return None, None, f"Failed to normalize PR2 mothur DB: {e}"

    if n_records <= 0:
        return None, None, "PR2 normalization produced zero records."
    return normalized_fasta, taxmap_tsv, None


def apply_pr2_taxmap_to_blast(blast_path, taxmap_path):
    if not os.path.exists(blast_path) or not os.path.exists(taxmap_path):
        return
    try:
        tax_df = pd.read_csv(taxmap_path, sep="\t", dtype=str)
        if tax_df.empty or "normalized_id" not in tax_df.columns:
            return
        tax_df = tax_df.fillna("")
        id_to_tax = dict(zip(tax_df["normalized_id"], tax_df.get("taxonomy", "")))
        id_to_orig = dict(zip(tax_df["normalized_id"], tax_df.get("original_id", "")))

        df = pd.read_csv(blast_path, sep="\t", header=None, names=BLAST_TABLE_COLUMNS, dtype=str)
        if df.empty:
            return

        def _mapped_title(sseqid, old_title):
            tax = str(id_to_tax.get(str(sseqid), "")).strip()
            orig = str(id_to_orig.get(str(sseqid), "")).strip()
            if tax and tax != "NA":
                return f"{tax} [{orig}]" if orig else tax
            if old_title and str(old_title).strip():
                return str(old_title)
            return orig or "NA"

        df["stitle"] = [
            _mapped_title(sseqid, old)
            for sseqid, old in zip(df["sseqid"], df["stitle"])
        ]
        df.to_csv(blast_path, sep="\t", header=False, index=False)
    except Exception:
        return

def ensure_blast_db(fasta_path, db_prefix):
    db_files = [f"{db_prefix}.nhr", f"{db_prefix}.nin", f"{db_prefix}.nsq"]
    if all(os.path.exists(p) for p in db_files):
        return True, ""
    cmd = [MAKEBLASTDB_PATH, "-in", fasta_path, "-dbtype", "nucl", "-parse_seqids", "-out", db_prefix]
    timeout_value = int(MAKEBLASTDB_TIMEOUT_SEC)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_value)
    except subprocess.TimeoutExpired:
        return False, f"makeblastdb timed out after {timeout_value}s for {os.path.basename(fasta_path)}"
    ok = all(os.path.exists(p) for p in db_files)
    return ok, result.stdout + result.stderr

def run_blastn(query_fasta, db_target, output_path, max_hits=5, evalue="1e-10", remote=False, timeout_sec=BLAST_TIMEOUT_SEC):
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle"
    cmd = [
        BLASTN_PATH,
        "-query", query_fasta,
        "-out", output_path,
        "-outfmt", outfmt,
        "-max_target_seqs", str(max_hits),
        "-evalue", str(evalue)
    ]
    if remote:
        cmd += ["-db", "nt", "-remote"]
    else:
        cmd += ["-db", db_target]
    timeout_value = int(timeout_sec) if timeout_sec is not None else int(BLAST_TIMEOUT_SEC)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_value)
    except subprocess.TimeoutExpired:
        return "", f"BLAST timed out after {timeout_value}s (db={db_target}, remote={remote}).", 124
    return result.stdout, result.stderr, result.returncode

def run_dada2(mode, output_dir, single_files=None, fwd_files=None, rev_files=None, ref_fasta=None, log_path=None, timeout_sec=None):
    script_path = os.path.join(os.getcwd(), "r", "run_dada2_pipeline.R")
    if not shutil.which(RSCRIPT_PATH):
        stderr = f"[ERROR] Rscript not found at {RSCRIPT_PATH}. Install R or set RSCRIPT_PATH."
        if log_path:
            _write_text_file(log_path, stderr)
        return "", stderr, 1
    if not os.path.exists(script_path):
        stderr = f"[ERROR] DADA2 pipeline script not found at {script_path}."
        if log_path:
            _write_text_file(log_path, stderr)
        return "", stderr, 1

    mode = str(mode or "").upper()
    if mode not in {"SE", "PE"}:
        stderr = f"[ERROR] Invalid DADA2 mode: {mode}."
        if log_path:
            _write_text_file(log_path, stderr)
        return "", stderr, 1

    if mode == "PE":
        if not fwd_files or not rev_files or len(fwd_files) != len(rev_files):
            stderr = "[ERROR] DADA2 PE mode requires aligned forward and reverse file lists."
            if log_path:
                _write_text_file(log_path, stderr)
            return "", stderr, 1
    if mode == "SE" and not single_files:
        stderr = "[ERROR] DADA2 SE mode requires single-end input files."
        if log_path:
            _write_text_file(log_path, stderr)
        return "", stderr, 1

    os.makedirs(output_dir, exist_ok=True)
    params_path = os.path.join(output_dir, "dada2_params.json")
    params_payload = {
        "mode": mode,
        "output_dir": output_dir,
        "single_files": list(single_files or []),
        "fwd_files": list(fwd_files or []),
        "rev_files": list(rev_files or []),
        "ref_fasta": ref_fasta or "",
    }
    _write_json(params_path, params_payload)

    cmd = [
        RSCRIPT_PATH,
        script_path,
        params_path,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec or BOLDIGGER_TIMEOUT_SEC)
    except subprocess.TimeoutExpired:
        stderr = f"DADA2 timed out after {timeout_sec}s" if timeout_sec else "DADA2 timed out."
        if log_path:
            _write_text_file(log_path, stderr)
        return "", stderr, 124

    stdout = result.stdout or ""
    stderr = result.stderr or ""
    if log_path:
        _write_text_file(
            log_path,
            "\n".join([
                "# DADA2 command",
                " ".join(cmd),
                "",
                "# stdout",
                stdout.strip(),
                "",
                "# stderr",
                stderr.strip(),
            ]) + "\n",
        )
    return stdout, stderr, result.returncode

def _safe_version(cmd):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return ""
        first_line = (result.stdout or result.stderr).splitlines()[0:1]
        return first_line[0] if first_line else ""
    except Exception:
        return ""

def _collect_tool_versions():
    return {
        "fastqc": _safe_version([FASTQC_PATH, "--version"]),
        "vsearch": _safe_version([VSEARCH_PATH, "--version"]),
        "mafft": _safe_version([MAFFT_PATH, "--version"]),
        "fasttree": _safe_version([FASTTREE_PATH, "-version"]),
        "raxml-ng": _safe_version([RAXML_PATH, "--version"]),
        "blastn": _safe_version([BLASTN_PATH, "-version"]),
        "Rscript": _safe_version([RSCRIPT_PATH, "--version"])
    }

def _compute_presence_jaccard_matrix(presence_df):
    if presence_df is None or presence_df.empty:
        return None

    matrix = presence_df.copy().fillna(0)
    samples = list(matrix.columns)
    if not samples:
        return None

    matrix = (matrix > 0).astype(int)
    values = {}
    for a in samples:
        values[a] = {}
        a_vec = matrix[a]
        for b in samples:
            b_vec = matrix[b]
            inter = int(((a_vec == 1) & (b_vec == 1)).sum())
            union = int(((a_vec == 1) | (b_vec == 1)).sum())
            values[a][b] = round(float(inter) / float(union), 4) if union else 0.0
    return pd.DataFrame(values).reindex(index=samples, columns=samples)


def _shared_and_unique_taxa_from_presence(presence_df):
    if presence_df is None or presence_df.empty:
        return [], {}

    binary = (presence_df.fillna(0) > 0).astype(int)
    n_samples = len(binary.columns)
    shared_taxa = binary.index[binary.sum(axis=1) == n_samples].tolist()

    unique_by_sample = {}
    for sample in binary.columns:
        unique_mask = (binary[sample] == 1) & (binary.drop(columns=[sample]).sum(axis=1) == 0)
        unique_by_sample[sample] = binary.index[unique_mask].tolist()

    return shared_taxa, unique_by_sample

def _select_blast_files(results_dir, db_filter=None):
    if not os.path.exists(results_dir):
        return []
    selected = None
    if db_filter:
        selected = set([str(x).strip().lower() for x in db_filter if str(x).strip()])
    blast_files = []
    for name in os.listdir(results_dir):
        if not (name.startswith("blast_") and name.endswith(".tsv")):
            continue
        db_name = name.replace("blast_", "").replace(".tsv", "").strip().lower()
        if selected and db_name not in selected:
            continue
        blast_files.append(os.path.join(results_dir, name))
    return sorted(blast_files)


def _build_blast_summary(results_dir, db_filter=None):
    if not os.path.exists(results_dir):
        return None
    blast_files = _select_blast_files(results_dir, db_filter=db_filter)
    if not blast_files:
        return None

    rows = []
    for blast_path in blast_files:
        db_name = os.path.basename(blast_path).replace("blast_", "").replace(".tsv", "")
        df = pd.read_csv(blast_path, sep="\t", header=None, names=BLAST_TABLE_COLUMNS)
        if df.empty:
            continue
        top = df.sort_values(["qseqid", "bitscore", "pident"], ascending=[True, False, False]).drop_duplicates("qseqid")
        top = top.copy()
        top["db"] = db_name
        rows.append(top[["db", "qseqid", "sseqid", "pident", "length", "qcovs", "evalue", "bitscore", "stitle"]])

    if not rows:
        return None

    return pd.concat(rows, ignore_index=True)


BLAST_TABLE_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "stitle"
]


def _clean_taxon_token(value):
    token = str(value or "").strip()
    if not token:
        return ""

    token = re.sub(r"^[a-zA-Z]__", "", token)
    token = re.sub(r"^[a-zA-Z]:", "", token)
    token = re.sub(r"\s+", " ", token).strip(" _")
    return token


def _normalize_taxon_label(taxon):
    raw = str(taxon or "").strip()
    if not raw:
        return ""

    parts = [_clean_taxon_token(p) for p in re.split(r"[;|]", raw)]
    parts = [p for p in parts if p]
    if not parts:
        return ""

    ignore_tokens = {
        "root", "cellular organisms", "eukaryota", "bacteria", "archaea",
        "unclassified", "uncultured", "unknown", "metagenome", "environmental sample", "na"
    }

    filtered = []
    for p in parts:
        lower = p.lower()
        if lower in ignore_tokens:
            continue
        if lower.startswith("uncultured") or lower.startswith("unclassified"):
            continue
        filtered.append(p)

    target = filtered[-1] if filtered else parts[-1]
    words = [w for w in target.split(" ") if w]
    if len(words) >= 2 and words[0][0].isalpha() and words[0][0].isupper():
        # Collapse binomial names to genus level for cross-DB stability.
        target = words[0]

    return target.lower()


def _extract_taxon_rank_labels(taxon):
    raw = str(taxon or "").strip()
    if not raw:
        return {"species": "", "genus": "", "family": ""}

    ignore_tokens = {
        "root", "cellular organisms", "eukaryota", "bacteria", "archaea",
        "unclassified", "uncultured", "unknown", "metagenome", "environmental sample", "na"
    }
    lineage = []
    for part in re.split(r"[;|]", raw):
        cleaned = _clean_taxon_token(part)
        if not cleaned:
            continue
        lower = cleaned.lower()
        if lower in ignore_tokens or lower.startswith("uncultured") or lower.startswith("unclassified"):
            continue
        lineage.append(cleaned)

    specific = lineage[-1] if lineage else raw
    words = [w for w in specific.split() if w]

    species = ""
    genus = ""
    if len(words) >= 2:
        genus = words[0].lower()
        species = f"{words[0].lower()} {words[1].lower()}"
    elif words:
        genus = words[0].lower()

    family = ""
    if len(lineage) >= 3:
        family = lineage[-3].lower()
    elif len(lineage) == 2:
        family = lineage[0].lower()
    elif genus:
        family = genus

    return {"species": species, "genus": genus, "family": family}


def _summarize_rank_agreement(rank_labels_by_db):
    available_levels = [level for level in ("species", "genus", "family") if [label for label in rank_labels_by_db.get(level, []) if label]]
    db_count = len(rank_labels_by_db.get("dbs", []))

    for level in ("species", "genus", "family"):
        labels = [label for label in rank_labels_by_db.get(level, []) if label]
        if len(labels) < 2:
            continue
        counts = pd.Series(labels).value_counts()
        top_label = str(counts.index[0])
        top_count = int(counts.iloc[0])
        if top_count >= 2:
            support = round(float(top_count) / float(len(labels)), 3) if labels else 0.0
            agreement = "full" if top_count == len(labels) and len(labels) >= 2 else "partial"
            return {
                "agreement_level": level,
                "agreement_support": support,
                "consensus_taxon_ranked": top_label,
                "agreement": agreement,
                "ranks_tested": ",".join(available_levels) if available_levels else "",
            }

    fallback_labels = [label for label in rank_labels_by_db.get("genus", []) if label] or [label for label in rank_labels_by_db.get("family", []) if label]
    counts = pd.Series(fallback_labels).value_counts() if fallback_labels else pd.Series(dtype=int)
    top_label = str(counts.index[0]) if not counts.empty else ""
    top_count = int(counts.iloc[0]) if not counts.empty else 0
    support = round(float(top_count) / float(len(fallback_labels)), 3) if fallback_labels else 0.0
    return {
        "agreement_level": "none",
        "agreement_support": support,
        "consensus_taxon_ranked": top_label,
        "agreement": "none",
        "ranks_tested": ",".join(available_levels) if available_levels else "",
    }


def _read_best_blast_hits(blast_path):
    try:
        df = pd.read_csv(blast_path, sep="\t", header=None, names=BLAST_TABLE_COLUMNS)
    except Exception:
        return pd.DataFrame(columns=["qseqid", "taxon", "bitscore", "pident", "evalue", "qcovs"])

    if df.empty:
        return pd.DataFrame(columns=["qseqid", "taxon", "bitscore", "pident", "evalue", "qcovs"])

    ranked = df.sort_values(["qseqid", "bitscore", "pident"], ascending=[True, False, False]).drop_duplicates("qseqid")
    out = ranked[["qseqid", "stitle", "bitscore", "pident", "evalue", "qcovs"]].copy()
    out["taxon"] = out["stitle"].astype(str).str.strip()
    out["taxon_norm"] = out["taxon"].apply(_normalize_taxon_label)
    rank_labels = out["taxon"].apply(_extract_taxon_rank_labels)
    out["taxon_species"] = rank_labels.apply(lambda value: value["species"])
    out["taxon_genus"] = rank_labels.apply(lambda value: value["genus"])
    out["taxon_family"] = rank_labels.apply(lambda value: value["family"])
    out = out.drop(columns=["stitle"])
    return out


def _build_blast_db_comparison(results_dir, db_filter=None):
    if not os.path.exists(results_dir):
        return None, None

    blast_files = _select_blast_files(results_dir, db_filter=db_filter)
    if len(blast_files) < 2:
        return None, None

    top_hits_by_db = {}
    db_hit_sets_exact = {}
    db_hit_sets_norm = {}
    all_qseqids = set()

    for blast_path in blast_files:
        db_name = os.path.basename(blast_path).replace("blast_", "").replace(".tsv", "")
        top_hits = _read_best_blast_hits(blast_path)
        if top_hits.empty:
            continue
        top_hits_by_db[db_name] = top_hits
        all_qseqids.update(top_hits["qseqid"].astype(str).tolist())
        db_hit_sets_exact[db_name] = set([
            (str(row["qseqid"]), str(row["taxon"])) for _, row in top_hits.iterrows()
        ])
        db_hit_sets_norm[db_name] = set([
            (str(row["qseqid"]), str(row["taxon_norm"]))
            for _, row in top_hits.iterrows()
            if str(row["taxon_norm"]).strip()
        ])

    if len(top_hits_by_db) < 2 or not all_qseqids:
        return None, None

    comparison_df = pd.DataFrame({"qseqid": sorted(all_qseqids)})
    for db_name, top_hits in top_hits_by_db.items():
        mapping = dict(zip(top_hits["qseqid"].astype(str), top_hits["taxon"].astype(str)))
        comparison_df[db_name] = comparison_df["qseqid"].map(mapping).fillna("")
        norm_mapping = dict(zip(top_hits["qseqid"].astype(str), top_hits["taxon_norm"].astype(str)))
        comparison_df[f"{db_name}_normalized"] = comparison_df["qseqid"].map(norm_mapping).fillna("")
        species_mapping = dict(zip(top_hits["qseqid"].astype(str), top_hits["taxon_species"].astype(str)))
        genus_mapping = dict(zip(top_hits["qseqid"].astype(str), top_hits["taxon_genus"].astype(str)))
        family_mapping = dict(zip(top_hits["qseqid"].astype(str), top_hits["taxon_family"].astype(str)))
        comparison_df[f"{db_name}_species"] = comparison_df["qseqid"].map(species_mapping).fillna("")
        comparison_df[f"{db_name}_genus"] = comparison_df["qseqid"].map(genus_mapping).fillna("")
        comparison_df[f"{db_name}_family"] = comparison_df["qseqid"].map(family_mapping).fillna("")

    db_columns = list(top_hits_by_db.keys())
    db_norm_columns = [f"{name}_normalized" for name in db_columns]
    support_values = []
    ratio_values = []
    unique_values = []
    consensus_taxa = []
    agreement_values = []
    db_with_hit_values = []
    agreement_level_values = []
    agreement_support_values = []
    ranks_tested_values = []
    consensus_ranked_values = []

    for _, row in comparison_df.iterrows():
        taxa = [str(row[c]).strip() for c in db_norm_columns if str(row[c]).strip()]
        if not taxa:
            support_values.append(0)
            ratio_values.append(0.0)
            unique_values.append(0)
            consensus_taxa.append("")
            agreement_values.append("none")
            db_with_hit_values.append(0)
            agreement_level_values.append("none")
            agreement_support_values.append(0.0)
            ranks_tested_values.append("")
            consensus_ranked_values.append("")
            continue

        counts = pd.Series(taxa).value_counts()
        support = int(counts.iloc[0])
        consensus = str(counts.index[0])
        db_with_hit = len(taxa)
        ratio = round(float(support) / float(db_with_hit), 3) if db_with_hit else 0.0
        unique_n = int(len(counts))

        db_rank_values = {"dbs": db_columns, "species": [], "genus": [], "family": []}
        for db_name in db_columns:
            if str(row.get(f"{db_name}_normalized", "")).strip():
                db_rank_values["species"].append(str(row.get(f"{db_name}_species", "")).strip())
                db_rank_values["genus"].append(str(row.get(f"{db_name}_genus", "")).strip())
                db_rank_values["family"].append(str(row.get(f"{db_name}_family", "")).strip())

        rank_summary = _summarize_rank_agreement(db_rank_values)
        agreement = rank_summary["agreement"]

        support_values.append(support)
        ratio_values.append(ratio)
        unique_values.append(unique_n)
        consensus_taxa.append(consensus)
        agreement_values.append(agreement)
        db_with_hit_values.append(db_with_hit)
        agreement_level_values.append(rank_summary["agreement_level"])
        agreement_support_values.append(rank_summary["agreement_support"])
        ranks_tested_values.append(rank_summary["ranks_tested"])
        consensus_ranked_values.append(rank_summary["consensus_taxon_ranked"])

    comparison_df["n_databases_with_hit"] = db_with_hit_values
    comparison_df["n_unique_taxa"] = unique_values
    comparison_df["consensus_taxon"] = consensus_taxa
    comparison_df["consensus_support"] = support_values
    comparison_df["consensus_ratio"] = ratio_values
    comparison_df["agreement"] = agreement_values
    comparison_df["agreement_level"] = agreement_level_values
    comparison_df["agreement_support"] = agreement_support_values
    comparison_df["ranks_tested"] = ranks_tested_values
    comparison_df["consensus_taxon_ranked"] = consensus_ranked_values

    pairwise_rows = []
    db_names = sorted(db_hit_sets_exact.keys())
    for i in range(len(db_names)):
        for j in range(i + 1, len(db_names)):
            a = db_names[i]
            b = db_names[j]
            set_a_exact = db_hit_sets_exact.get(a, set())
            set_b_exact = db_hit_sets_exact.get(b, set())
            inter_exact = len(set_a_exact.intersection(set_b_exact))
            union_exact = len(set_a_exact.union(set_b_exact))
            jaccard_exact = round(float(inter_exact) / float(union_exact), 4) if union_exact else 0.0

            set_a_norm = db_hit_sets_norm.get(a, set())
            set_b_norm = db_hit_sets_norm.get(b, set())
            inter_norm = len(set_a_norm.intersection(set_b_norm))
            union_norm = len(set_a_norm.union(set_b_norm))
            jaccard_norm = round(float(inter_norm) / float(union_norm), 4) if union_norm else 0.0

            pairwise_rows.append({
                "db_a": a,
                "db_b": b,
                "shared_assignments_exact": inter_exact,
                "union_assignments_exact": union_exact,
                "jaccard_exact": jaccard_exact,
                "shared_assignments_normalized": inter_norm,
                "union_assignments_normalized": union_norm,
                "jaccard_normalized": jaccard_norm
            })

    pairwise_df = pd.DataFrame(pairwise_rows) if pairwise_rows else None
    return comparison_df, pairwise_df


def _to_numeric(series):
    return pd.to_numeric(series, errors="coerce")


def _build_db_metrics_table(results_dir, comparison_df=None, pairwise_df=None, db_filter=None):
    blast_files = _select_blast_files(results_dir, db_filter=db_filter)
    if not blast_files:
        return None

    rows = []
    for blast_path in blast_files:
        db_name = os.path.basename(blast_path).replace("blast_", "").replace(".tsv", "")
        best_hits = _read_best_blast_hits(blast_path)
        if best_hits.empty:
            rows.append({
                "db": db_name,
                "n_queries_with_hit": 0,
                "n_unique_taxa_normalized": 0,
                "mean_identity_pct": 0.0,
                "median_bitscore": 0.0,
                "median_evalue": None,
                "mean_consensus_ratio": 0.0,
                "full_agreement_queries": 0,
                "partial_agreement_queries": 0,
                "none_agreement_queries": 0,
                "avg_pairwise_jaccard_normalized": 0.0,
            })
            continue

        n_queries = int(best_hits["qseqid"].astype(str).nunique())
        n_taxa_norm = int(best_hits["taxon_norm"].astype(str).str.strip().replace("", pd.NA).dropna().nunique())
        pident = _to_numeric(best_hits["pident"])
        bitscore = _to_numeric(best_hits["bitscore"])
        evalue = _to_numeric(best_hits["evalue"])

        mean_identity = round(float(pident.mean()), 3) if not pident.dropna().empty else 0.0
        median_bitscore = round(float(bitscore.median()), 3) if not bitscore.dropna().empty else 0.0
        median_evalue = float(evalue.median()) if not evalue.dropna().empty else None

        mean_consensus_ratio = 0.0
        full_count = 0
        partial_count = 0
        none_count = 0
        db_norm_col = f"{db_name}_normalized"
        if comparison_df is not None and not comparison_df.empty and db_norm_col in comparison_df.columns:
            subset = comparison_df[comparison_df[db_norm_col].astype(str).str.strip() != ""]
            if not subset.empty:
                ratio_values = _to_numeric(subset.get("consensus_ratio", pd.Series(dtype=float)))
                mean_consensus_ratio = round(float(ratio_values.mean()), 4) if not ratio_values.dropna().empty else 0.0
                agreements = subset.get("agreement", pd.Series(dtype=str)).astype(str).str.strip().str.lower()
                full_count = int((agreements == "full").sum())
                partial_count = int((agreements == "partial").sum())
                none_count = int((agreements == "none").sum())

        avg_pairwise = 0.0
        if pairwise_df is not None and not pairwise_df.empty:
            mask = (pairwise_df["db_a"] == db_name) | (pairwise_df["db_b"] == db_name)
            subset_pairs = pairwise_df[mask]
            if not subset_pairs.empty:
                jvals = _to_numeric(subset_pairs.get("jaccard_normalized", pd.Series(dtype=float)))
                avg_pairwise = round(float(jvals.mean()), 4) if not jvals.dropna().empty else 0.0

        rows.append({
            "db": db_name,
            "n_queries_with_hit": n_queries,
            "n_unique_taxa_normalized": n_taxa_norm,
            "mean_identity_pct": mean_identity,
            "median_bitscore": median_bitscore,
            "median_evalue": median_evalue,
            "mean_consensus_ratio": mean_consensus_ratio,
            "full_agreement_queries": full_count,
            "partial_agreement_queries": partial_count,
            "none_agreement_queries": none_count,
            "avg_pairwise_jaccard_normalized": avg_pairwise,
        })

    metrics_df = pd.DataFrame(rows)
    if metrics_df.empty:
        return None
    return metrics_df.sort_values(["mean_consensus_ratio", "avg_pairwise_jaccard_normalized", "n_queries_with_hit"], ascending=[False, False, False]).reset_index(drop=True)


def _build_analysis_overview_table(presence_df, db_metrics_df):
    rows = []

    if db_metrics_df is not None and not db_metrics_df.empty:
        n_dbs = int(len(db_metrics_df))
        mean_hits = round(float(_to_numeric(db_metrics_df["n_queries_with_hit"]).mean()), 2)
        mean_consensus = round(float(_to_numeric(db_metrics_df["mean_consensus_ratio"]).mean()), 4)
        mean_pairwise = round(float(_to_numeric(db_metrics_df["avg_pairwise_jaccard_normalized"]).mean()), 4)
        rows.append({
            "analysis": "blast_multi_db",
            "n_entities": n_dbs,
            "n_features": int(_to_numeric(db_metrics_df["n_unique_taxa_normalized"]).sum()),
            "shared_features": None,
            "mean_unique_features": None,
            "quality_indicator": f"mean_hits={mean_hits}; mean_consensus={mean_consensus}; mean_pairwise_jaccard={mean_pairwise}",
        })

    if not rows:
        return None
    return pd.DataFrame(rows)


def _load_otu_abundance_matrix(otu_matrix_path):
    if not otu_matrix_path or not os.path.exists(otu_matrix_path) or os.path.getsize(otu_matrix_path) == 0:
        return None
    try:
        matrix = pd.read_csv(otu_matrix_path, sep="\t", index_col=0)
    except Exception:
        return None
    if matrix.empty:
        return None
    return matrix.apply(pd.to_numeric, errors="coerce").fillna(0)


def _calculate_alpha_diversity(otu_matrix):
    if otu_matrix is None or otu_matrix.empty:
        return None

    rows = []
    for sample in otu_matrix.columns:
        counts = pd.to_numeric(otu_matrix[sample], errors="coerce").fillna(0).astype(float)
        total = float(counts.sum())
        observed = int((counts > 0).sum())
        if total > 0:
            proportions = counts[counts > 0] / total
            shannon = float(-(proportions * np.log(proportions)).sum()) if not proportions.empty else 0.0
            simpson = float(1.0 - (proportions ** 2).sum()) if not proportions.empty else 0.0
        else:
            shannon = 0.0
            simpson = 0.0

        singletons = int((counts == 1).sum())
        doubletons = int((counts == 2).sum())
        if doubletons > 0:
            chao1 = observed + ((singletons ** 2) / (2.0 * doubletons))
        else:
            chao1 = observed + (singletons * (singletons - 1) / 2.0)

        rows.append({
            "sample": sample,
            "richness": observed,
            "shannon": round(shannon, 6),
            "simpson": round(simpson, 6),
            "chao1": round(float(chao1), 6),
            "total_abundance": int(total),
        })

    return pd.DataFrame(rows).set_index("sample")


def _pcoa_from_distance(distance_df):
    if distance_df is None or distance_df.empty:
        return None

    distances = distance_df.to_numpy(dtype=float)
    n = distances.shape[0]
    if n < 2:
        return None

    h = np.eye(n) - np.ones((n, n)) / float(n)
    b = -0.5 * h @ (distances ** 2) @ h
    eigenvalues, eigenvectors = np.linalg.eigh(b)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    positive = eigenvalues > 0
    if not positive.any():
        return None

    dims = min(2, int(positive.sum()))
    coords = eigenvectors[:, positive][:, :dims] * np.sqrt(eigenvalues[positive][:dims])
    columns = [f"PCoA{i}" for i in range(1, dims + 1)]
    return pd.DataFrame(coords, index=distance_df.index, columns=columns)


def _calculate_beta_diversity(otu_matrix):
    if otu_matrix is None or otu_matrix.empty:
        return None, None, None, None

    sample_names = list(otu_matrix.columns)
    presence = (otu_matrix > 0).astype(int)
    jaccard = pd.DataFrame(index=sample_names, columns=sample_names, dtype=float)
    bray_curtis = pd.DataFrame(index=sample_names, columns=sample_names, dtype=float)

    for a in sample_names:
        a_presence = presence[a].to_numpy(dtype=int)
        a_counts = otu_matrix[a].to_numpy(dtype=float)
        for b in sample_names:
            b_presence = presence[b].to_numpy(dtype=int)
            b_counts = otu_matrix[b].to_numpy(dtype=float)
            intersection = float(np.logical_and(a_presence > 0, b_presence > 0).sum())
            union = float(np.logical_or(a_presence > 0, b_presence > 0).sum())
            jaccard.at[a, b] = 1.0 - (intersection / union) if union else 0.0

            numerator = float(np.abs(a_counts - b_counts).sum())
            denominator = float((a_counts + b_counts).sum())
            bray_curtis.at[a, b] = numerator / denominator if denominator else 0.0

    pcoa_coords = _pcoa_from_distance(jaccard)

    nmds_coords = None
    try:
        from sklearn.manifold import MDS
        nmds_model = MDS(n_components=2, metric=False, dissimilarity="precomputed", random_state=42, n_init=4, max_iter=300)
        nmds_array = nmds_model.fit_transform(jaccard.to_numpy(dtype=float))
        nmds_coords = pd.DataFrame(nmds_array, index=sample_names, columns=["NMDS1", "NMDS2"])
    except Exception:
        nmds_coords = None

    return jaccard, bray_curtis, pcoa_coords, nmds_coords


def _build_shared_unique_taxa(otu_matrix):
    if otu_matrix is None or otu_matrix.empty:
        return None, None

    presence = (otu_matrix > 0).astype(int)
    shared_taxa, unique_by_sample = _shared_and_unique_taxa_from_presence(presence)

    shared_df = pd.DataFrame({"otu_id": shared_taxa})
    unique_rows = []
    for sample, taxa in unique_by_sample.items():
        for taxon in taxa:
            unique_rows.append({"sample": sample, "otu_id": taxon})
    unique_df = pd.DataFrame(unique_rows)
    return shared_df, unique_df


def _build_heatmap_abundance(otu_matrix):
    if otu_matrix is None or otu_matrix.empty:
        return None
    normalized = otu_matrix.div(otu_matrix.sum(axis=0).replace(0, np.nan), axis=1).fillna(0) * 1_000_000.0
    return np.log1p(normalized)


def _write_dataframe_tsv(df, path, index=True, index_label=None):
    if df is None or path is None:
        return None
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        df.to_csv(path, sep="\t", index=index, index_label=index_label)
    except Exception:
        return None
    return path if _path_is_valid(path) else None


def _write_plotly_html(fig, path):
    if fig is None or path is None:
        return None
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fig.write_html(path, include_plotlyjs="cdn", full_html=True)
    except Exception:
        return None
    return path if _path_is_valid(path) else None


def _build_stacked_abundance_table(otu_matrix, top_n=20):
    if otu_matrix is None or otu_matrix.empty:
        return None, None

    matrix = otu_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
    totals = matrix.sum(axis=1).sort_values(ascending=False)
    keep_features = list(totals.head(max(1, int(top_n))).index)
    selected = matrix.loc[keep_features].copy()
    other = matrix.drop(index=keep_features, errors="ignore")
    if not other.empty:
        selected.loc["Other"] = other.sum(axis=0)

    if selected.empty:
        return None, None

    long_rows = []
    for sample in selected.columns:
        sample_total = float(selected[sample].sum())
        for feature in selected.index:
            abundance = float(selected.at[feature, sample])
            long_rows.append({
                "sample": sample,
                "feature": feature,
                "abundance": abundance,
                "relative_abundance": (abundance / sample_total) if sample_total else 0.0,
                "sample_total": sample_total,
            })

    long_df = pd.DataFrame(long_rows)
    fig = go.Figure()
    for feature in selected.index:
        fig.add_trace(go.Bar(
            name=str(feature),
            x=list(selected.columns),
            y=selected.loc[feature].astype(float).tolist(),
        ))
    fig.update_layout(
        barmode="stack",
        title="Stacked abundance by sample",
        xaxis_title="Sample",
        yaxis_title="Reads",
        legend_title="Feature",
        template="plotly_white",
    )
    return long_df, fig


def _build_rarefaction_curves(otu_matrix, max_points=12, replicates=5):
    if otu_matrix is None or otu_matrix.empty:
        return None, None

    matrix = otu_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
    rows = []
    fig = go.Figure()

    for sample in matrix.columns:
        counts = matrix[sample].astype(int).clip(lower=0)
        feature_labels = []
        for feature, count in counts.items():
            if count > 0:
                feature_labels.extend([str(feature)] * int(count))
        total_reads = len(feature_labels)
        if total_reads == 0:
            continue

        depths = np.unique(np.linspace(1, total_reads, num=min(int(max_points), total_reads), dtype=int))
        curve_replicates = []
        seed = abs(hash(sample)) % (2**32)
        rng = np.random.default_rng(seed)

        for rep in range(max(1, int(replicates))):
            permuted = np.array(feature_labels, dtype=object)
            rng.shuffle(permuted)
            seen = set()
            richness_values = []
            cursor = 0
            for depth in depths:
                while cursor < depth and cursor < total_reads:
                    seen.add(permuted[cursor])
                    cursor += 1
                richness_values.append(float(len(seen)))
            curve_replicates.append(richness_values)

        curve_array = np.array(curve_replicates, dtype=float)
        mean_curve = curve_array.mean(axis=0).tolist()
        std_curve = curve_array.std(axis=0).tolist() if len(curve_replicates) > 1 else [0.0] * len(mean_curve)

        for depth, richness_mean, richness_std in zip(depths, mean_curve, std_curve):
            rows.append({
                "sample": sample,
                "depth": int(depth),
                "richness_mean": round(float(richness_mean), 6),
                "richness_std": round(float(richness_std), 6),
                "total_reads": int(total_reads),
            })

        fig.add_trace(go.Scatter(
            x=depths.tolist(),
            y=mean_curve,
            mode="lines+markers",
            name=str(sample),
        ))

    if not rows:
        return None, None

    fig.update_layout(
        title="Rarefaction curves",
        xaxis_title="Sequencing depth",
        yaxis_title="Observed richness",
        template="plotly_white",
    )
    return pd.DataFrame(rows), fig


def _export_analysis_outputs(project_state, method_name, otu_matrix):
    if otu_matrix is None or otu_matrix.empty:
        return {}

    method_slug = "otu" if str(method_name).lower() == "vsearch" else "asv" if str(method_name).lower() == "dada2" else _safe_name(method_name or "analysis")
    outputs = {}
    canonical_root = _step_outputs_base(project_state, "analysis")
    results_root = os.path.join(_get_active_root(project_state), "results")

    def _canon_path(rel_path):
        return os.path.join(canonical_root, rel_path)

    def _alias_path(rel_path):
        return os.path.join(results_root, rel_path)

    matrix_rel = f"{method_slug}_sample_matrix.tsv"
    matrix_path = _write_dataframe_tsv(otu_matrix, _canon_path(matrix_rel), index=True, index_label=method_slug.upper())
    if matrix_path:
        outputs[f"{method_slug}_sample_matrix"] = matrix_path
        outputs[f"{method_slug}_sample_matrix_result"] = _mirror_output_to_results(matrix_path, _alias_path(matrix_rel))

    alpha_df = _calculate_alpha_diversity(otu_matrix)
    if alpha_df is not None and not alpha_df.empty:
        alpha_rel = "alpha_diversity.tsv"
        alpha_path = _write_dataframe_tsv(alpha_df.reset_index(), _canon_path(alpha_rel), index=False)
        if alpha_path:
            outputs["alpha_diversity"] = alpha_path
            outputs["alpha_diversity_result"] = _mirror_output_to_results(alpha_path, _alias_path(alpha_rel))

    jaccard_df, bray_df, pcoa_df, nmds_df = _calculate_beta_diversity(otu_matrix)
    beta_rel = "beta_diversity"
    if jaccard_df is not None and not jaccard_df.empty:
        jaccard_path = _write_dataframe_tsv(jaccard_df, _canon_path(os.path.join(beta_rel, "jaccard_distance.tsv")))
        if jaccard_path:
            outputs["jaccard_distance"] = jaccard_path
            outputs["jaccard_distance_result"] = _mirror_output_to_results(jaccard_path, _alias_path(os.path.join(beta_rel, "jaccard_distance.tsv")))
    if bray_df is not None and not bray_df.empty:
        bray_path = _write_dataframe_tsv(bray_df, _canon_path(os.path.join(beta_rel, "braycurtis_distance.tsv")))
        if bray_path:
            outputs["braycurtis_distance"] = bray_path
            outputs["braycurtis_distance_result"] = _mirror_output_to_results(bray_path, _alias_path(os.path.join(beta_rel, "braycurtis_distance.tsv")))
    if pcoa_df is not None and not pcoa_df.empty:
        pcoa_path = _write_dataframe_tsv(pcoa_df, _canon_path(os.path.join(beta_rel, "pcoa_coords.tsv")))
        if pcoa_path:
            outputs["pcoa_coords"] = pcoa_path
            outputs["pcoa_coords_result"] = _mirror_output_to_results(pcoa_path, _alias_path(os.path.join(beta_rel, "pcoa_coords.tsv")))
    if nmds_df is not None and not nmds_df.empty:
        nmds_path = _write_dataframe_tsv(nmds_df, _canon_path(os.path.join(beta_rel, "nmds_coords.tsv")))
        if nmds_path:
            outputs["nmds_coords"] = nmds_path
            outputs["nmds_coords_result"] = _mirror_output_to_results(nmds_path, _alias_path(os.path.join(beta_rel, "nmds_coords.tsv")))

    shared_df, unique_df = _build_shared_unique_taxa(otu_matrix)
    comparative_rel = "comparative_analysis"
    if shared_df is not None and not shared_df.empty:
        shared_path = _write_dataframe_tsv(shared_df, _canon_path(os.path.join(comparative_rel, "shared_taxa.tsv")), index=False)
        if shared_path:
            outputs["shared_taxa"] = shared_path
            outputs["shared_taxa_result"] = _mirror_output_to_results(shared_path, _alias_path(os.path.join(comparative_rel, "shared_taxa.tsv")))
    if unique_df is not None and not unique_df.empty:
        unique_path = _write_dataframe_tsv(unique_df, _canon_path(os.path.join(comparative_rel, "unique_taxa_per_sample.tsv")), index=False)
        if unique_path:
            outputs["unique_taxa_per_sample"] = unique_path
            outputs["unique_taxa_per_sample_result"] = _mirror_output_to_results(unique_path, _alias_path(os.path.join(comparative_rel, "unique_taxa_per_sample.tsv")))

    heatmap_df = _build_heatmap_abundance(otu_matrix)
    if heatmap_df is not None and not heatmap_df.empty:
        heatmap_rel = os.path.join(comparative_rel, "heatmap_abundance.tsv")
        heatmap_path = _write_dataframe_tsv(heatmap_df, _canon_path(heatmap_rel))
        if heatmap_path:
            outputs["heatmap_abundance"] = heatmap_path
            outputs["heatmap_abundance_result"] = _mirror_output_to_results(heatmap_path, _alias_path(heatmap_rel))
        fig = go.Figure(data=go.Heatmap(
            z=heatmap_df.to_numpy(dtype=float),
            x=[str(c) for c in heatmap_df.columns],
            y=[str(i) for i in heatmap_df.index],
            colorscale="Viridis",
        ))
        fig.update_layout(title="Heatmap abundance", xaxis_title="Sample", yaxis_title="Feature", template="plotly_white")
        heatmap_html = _write_plotly_html(fig, _canon_path(os.path.join(comparative_rel, "heatmap_abundance.html")))
        if heatmap_html:
            outputs["heatmap_abundance_html"] = heatmap_html
            outputs["heatmap_abundance_html_result"] = _mirror_output_to_results(heatmap_html, _alias_path(os.path.join(comparative_rel, "heatmap_abundance.html")))

    stacked_df, stacked_fig = _build_stacked_abundance_table(otu_matrix)
    if stacked_df is not None and not stacked_df.empty:
        stacked_rel = os.path.join(comparative_rel, "stacked_abundance.tsv")
        stacked_path = _write_dataframe_tsv(stacked_df, _canon_path(stacked_rel), index=False)
        if stacked_path:
            outputs["stacked_abundance"] = stacked_path
            outputs["stacked_abundance_result"] = _mirror_output_to_results(stacked_path, _alias_path(stacked_rel))
        stacked_html = _write_plotly_html(stacked_fig, _canon_path(os.path.join(comparative_rel, "stacked_abundance.html")))
        if stacked_html:
            outputs["stacked_abundance_html"] = stacked_html
            outputs["stacked_abundance_html_result"] = _mirror_output_to_results(stacked_html, _alias_path(os.path.join(comparative_rel, "stacked_abundance.html")))

    rarefaction_df, rarefaction_fig = _build_rarefaction_curves(otu_matrix)
    if rarefaction_df is not None and not rarefaction_df.empty:
        rarefaction_rel = os.path.join(comparative_rel, "rarefaction_curves.tsv")
        rarefaction_path = _write_dataframe_tsv(rarefaction_df, _canon_path(rarefaction_rel), index=False)
        if rarefaction_path:
            outputs["rarefaction_curves"] = rarefaction_path
            outputs["rarefaction_curves_result"] = _mirror_output_to_results(rarefaction_path, _alias_path(rarefaction_rel))
        rarefaction_html = _write_plotly_html(rarefaction_fig, _canon_path(os.path.join(comparative_rel, "rarefaction_curves.html")))
        if rarefaction_html:
            outputs["rarefaction_curves_html"] = rarefaction_html
            outputs["rarefaction_curves_html_result"] = _mirror_output_to_results(rarefaction_html, _alias_path(os.path.join(comparative_rel, "rarefaction_curves.html")))

    if alpha_df is not None and not alpha_df.empty:
        alpha_reset = alpha_df.reset_index().rename(columns={"index": "sample"})
        sample_summary = alpha_reset.copy()
    else:
        sample_summary = pd.DataFrame(columns=["sample", "richness", "shannon", "simpson", "chao1", "total_abundance"])

    total_reads = int(float(otu_matrix.to_numpy(dtype=float).sum())) if otu_matrix is not None else 0
    run_summary_row = pd.DataFrame([{
        "scope": "run",
        "sample": "analysis",
        "n_samples": int(len(otu_matrix.columns)),
        "n_features": int(len(otu_matrix.index)),
        "total_reads": total_reads,
        "mean_reads_per_sample": round(float(total_reads / max(1, len(otu_matrix.columns))), 6),
        "matrix_kind": method_slug,
    }])
    if not sample_summary.empty:
        sample_summary = sample_summary.copy()
        sample_summary.insert(0, "scope", "sample")
        summary_df = pd.concat([run_summary_row, sample_summary], ignore_index=True, sort=False)
    else:
        summary_df = run_summary_row

    summary_rel = "summary_stats.tsv"
    summary_path = _write_dataframe_tsv(summary_df, _canon_path(summary_rel), index=False)
    if summary_path:
        outputs["summary_stats"] = summary_path
        outputs["summary_stats_result"] = _mirror_output_to_results(summary_path, _alias_path(summary_rel))

    return outputs

def fastq_to_fasta(fastq_path, fasta_path):
    if not os.path.exists(fastq_path) or os.path.getsize(fastq_path) == 0:
        return False
    opener = gzip.open if str(fastq_path).lower().endswith(".gz") else open
    with opener(fastq_path, "rt", encoding="utf-8", errors="replace") as fq, open(fasta_path, "w", encoding="utf-8") as fa:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline()
            plus = fq.readline()
            qual = fq.readline()
            if not (header and seq and plus and qual):
                break
            if header.startswith('@'):
                fa.write(">" + header[1:])
                fa.write(seq)
    return os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 0

def run_vsearch_otus(input_fasta, output_dir, vsearch_path=VSEARCH_PATH, identity=0.97):
    # Delegate to vsearch_pipeline implementation which returns a dict
    res = vp.run_vsearch_otus(input_fasta, output_dir, vsearch_path=vsearch_path, identity=identity)
    if not res or res.get('error'):
        logs = res.get('logs') if isinstance(res, dict) else {}
        derep_out = logs.get('derep','')
        cluster_out = logs.get('cluster','')
        return derep_out, cluster_out, res.get('repseqs') if res else '', res.get('uc_file') if res else ''
    # return outputs compatible with previous callers
    logs = res.get('logs', {})
    derep_out = logs.get('derep','')
    cluster_out = logs.get('cluster','')
    return derep_out, cluster_out, res.get('repseqs'), res.get('uc_file')


def _analysis_worker(analysis_mode, identity, uploaded_files, project_state, samples_store, workflow_scope, run_dir):
    """Background worker that executes the analysis pipeline."""
    try:
        ps = dict(project_state or {})
        ps['run_dir'] = run_dir

        methods = _analysis_mode_to_methods(analysis_mode)
        manifest = {"methods": {}, "created_at": datetime.now().isoformat(), "status": STEP_STATUS_RUNNING, "errors": []}
        analysis_artifacts = {}
        compat_results_root = os.path.join(_get_active_root(ps), "results")

        def _record_method_status(method_name, status, outputs=None, inputs=None, error=None, metrics=None):
            manifest["methods"][method_name] = {
                "status": status,
                "outputs": outputs or {},
                "inputs": inputs or {},
            }
            if error:
                manifest["methods"][method_name]["error"] = error
                manifest.setdefault("errors", []).append(error)
            if metrics is not None:
                manifest["methods"][method_name]["metrics"] = metrics

        for m in methods:
            if m == "vsearch":
                pooled_input_dir = _resolve_output_path(ps, "analysis", "intermediate", "pooled_analysis")
                fasta_file, sample_counts = _build_pooled_fasta_for_analysis(samples_store, workflow_scope, pooled_input_dir)
                if not fasta_file or not os.path.exists(fasta_file) or os.path.getsize(fasta_file) == 0:
                    _record_method_status("vsearch", STEP_STATUS_FAILED, outputs={}, inputs={"sample_ids": list((workflow_scope or {}).get("samples_selected", []))}, error="No valid pooled FASTA could be built for VSEARCH.")
                    continue
                if not sample_counts:
                    _record_method_status("vsearch", STEP_STATUS_FAILED, outputs={}, inputs={"sample_ids": list((workflow_scope or {}).get("samples_selected", []))}, error="No readable sample FASTA inputs for VSEARCH.")
                    continue

                output_dir_vsearch = os.path.join(_step_outputs_base(ps, "analysis"), "vsearch")
                os.makedirs(output_dir_vsearch, exist_ok=True)
                res = vp.run_vsearch_otus(fasta_file, output_dir_vsearch, vsearch_path=VSEARCH_PATH, identity=identity)
                if not res or res.get('error'):
                    _record_method_status("vsearch", STEP_STATUS_FAILED, outputs={}, inputs={"sample_ids": list((workflow_scope or {}).get("samples_selected", []))}, error=res.get('error') if isinstance(res, dict) else "vsearch_failed")
                    continue

                otus_fasta = res.get('repseqs')
                uc_file = res.get('uc_file')
                otu_representatives_path = _resolve_output_path(ps, "analysis", "final", "otus.fasta")
                try:
                    shutil.copy2(otus_fasta, otu_representatives_path)
                except Exception:
                    otu_representatives_path = otus_fasta
                _mirror_output_to_results(otu_representatives_path, os.path.join(compat_results_root, "otus.fasta"))

                pooled_matrix = _parse_vsearch_uc_to_matrix(uc_file, sample_counts) if uc_file else None
                if pooled_matrix is None or getattr(pooled_matrix, 'empty', False):
                    _record_method_status("vsearch", STEP_STATUS_FAILED, outputs={"otus_fasta": otu_representatives_path, "uc_file": uc_file}, inputs={"sample_ids": list((workflow_scope or {}).get("samples_selected", []))}, error="No per-sample OTU matrix could be built.")
                    continue

                vmetrics = {'n_otus': _count_fasta_seqs(otus_fasta) if os.path.exists(otus_fasta) else 0, 'samples': len(pooled_matrix.columns)}
                analysis_artifacts.update(_export_analysis_outputs(ps, "vsearch", pooled_matrix))
                _record_method_status("vsearch", STEP_STATUS_SUCCESS, outputs={"otus_fasta": otu_representatives_path, "uc_file": uc_file, 'otu_sample_matrix': analysis_artifacts.get('otu_sample_matrix')}, inputs={"sample_ids": list((workflow_scope or {}).get("samples_selected", []))}, metrics=vmetrics)

            elif m == "dada2":
                dada2_dir = os.path.join(_step_outputs_base(ps, "analysis"), "dada2")
                os.makedirs(dada2_dir, exist_ok=True)
                dada2_logs_dir = os.path.join(dada2_dir, "logs")
                os.makedirs(dada2_logs_dir, exist_ok=True)
                seqtab_out = os.path.join(dada2_dir, "seqtab_nochim.rds")
                asv_fasta_out = os.path.join(dada2_dir, "repseqs.fasta")
                seqtab_tsv_out = os.path.join(dada2_dir, "seqtab_nochim.csv")
                dada2_log_path = os.path.join(dada2_logs_dir, "dada2.log")

                dada2_inputs = _resolve_dada2_inputs(ps, uploaded_files, samples_store, workflow_scope)
                input_mode = dada2_inputs.get("mode")
                input_source = dada2_inputs.get("source")
                input_error = dada2_inputs.get("error")
                fwd_files = dada2_inputs.get("fwd_files") or []
                rev_files = dada2_inputs.get("rev_files") or []
                single_files = dada2_inputs.get("single_files") or []

                if input_error:
                    _record_method_status("dada2", STEP_STATUS_FAILED, outputs={}, inputs={"source": input_source or "none", "mode": input_mode or "unknown"}, error=input_error)
                    continue

                stdout, stderr, code = run_dada2(
                    input_mode,
                    dada2_dir,
                    single_files=single_files,
                    fwd_files=fwd_files,
                    rev_files=rev_files,
                    ref_fasta=DEFAULT_PR2_DADA2_FASTA_PATH if os.path.exists(DEFAULT_PR2_DADA2_FASTA_PATH) else "",
                    log_path=dada2_log_path,
                    timeout_sec=BOLDIGGER_TIMEOUT_SEC,
                )

                if code != 0 or not os.path.exists(asv_fasta_out) or not os.path.exists(seqtab_out):
                    _record_method_status("dada2", STEP_STATUS_FAILED, outputs={"asv_fasta": asv_fasta_out, "seqtab_nochim_rds": seqtab_out, 'log': dada2_log_path}, inputs={'mode': input_mode, 'source': input_source, 'fwd_files': fwd_files, 'rev_files': rev_files, 'single_files': single_files}, error="DADA2 did not produce expected outputs.")
                    continue

                final_asv = _resolve_output_path(ps, "analysis", "final", "asv.fasta")
                final_seqtab = _resolve_output_path(ps, "analysis", "final", "seqtab_nochim.rds")
                _mirror_output_to_results(asv_fasta_out, final_asv)
                _mirror_output_to_results(seqtab_out, final_seqtab)
                _mirror_output_to_results(final_asv, os.path.join(compat_results_root, "asv.fasta"))
                _mirror_output_to_results(final_seqtab, os.path.join(compat_results_root, "seqtab_nochim.rds"))

                seqtab_matrix = None
                if os.path.exists(seqtab_tsv_out) and os.path.getsize(seqtab_tsv_out) > 0:
                    try:
                        seqtab_matrix = pd.read_csv(seqtab_tsv_out)
                        if "OTU" in seqtab_matrix.columns:
                            seqtab_matrix = seqtab_matrix.set_index("OTU")
                        elif seqtab_matrix.columns.size > 0:
                            seqtab_matrix = seqtab_matrix.set_index(seqtab_matrix.columns[0])
                        seqtab_matrix = seqtab_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
                    except Exception:
                        seqtab_matrix = None
                if seqtab_matrix is None or getattr(seqtab_matrix, 'empty', False):
                    manifest.setdefault("warnings", []).append("DADA2 TSV output was not available or could not be parsed; downstream analysis exports were skipped.")

                metrics_path = os.path.join(dada2_dir, 'dada2_metrics.json')
                dmetrics = {}
                if os.path.exists(metrics_path):
                    try:
                        dmetrics = _read_json(metrics_path) or {}
                    except Exception:
                        dmetrics = {}
                else:
                    dmetrics['n_asvs'] = _count_fasta_seqs(asv_fasta_out) if os.path.exists(asv_fasta_out) else 0

                if seqtab_matrix is not None and not seqtab_matrix.empty:
                    analysis_artifacts.update(_export_analysis_outputs(ps, "dada2", seqtab_matrix))

                _record_method_status("dada2", STEP_STATUS_SUCCESS, outputs={'asv_fasta': asv_fasta_out, 'seqtab_nochim_rds': seqtab_out, 'seqtab_nochim_tsv': seqtab_tsv_out, 'asv_fasta_result': final_asv, 'seqtab_nochim_rds_result': final_seqtab, 'log': dada2_log_path}, inputs={'mode': input_mode, 'source': input_source, 'fwd_files': fwd_files, 'rev_files': rev_files, 'single_files': single_files}, metrics=dmetrics)

        if analysis_mode == "both":
            v_ok = (manifest.get("methods", {}).get("vsearch", {}) or {}).get("status") == STEP_STATUS_SUCCESS
            d_ok = (manifest.get("methods", {}).get("dada2", {}) or {}).get("status") == STEP_STATUS_SUCCESS
            if v_ok and d_ok:
                comparison_dir = os.path.join(_step_outputs_base(ps, "analysis"), "comparison")
                os.makedirs(comparison_dir, exist_ok=True)
                summary_payload = {
                    "mode": "both",
                    "created_at": datetime.now().isoformat(timespec="seconds"),
                    "vsearch": manifest.get("methods", {}).get("vsearch", {}).get("metrics", {}),
                    "dada2": manifest.get("methods", {}).get("dada2", {}).get("metrics", {}),
                    "note": "ASV and OTU outputs are kept independent; comparison is derived from metrics only.",
                }
                summary_path = os.path.join(comparison_dir, "summary.json")
                _write_json(summary_path, summary_payload)
                analysis_artifacts["comparison_summary"] = summary_path
                analysis_artifacts["comparison_summary_result"] = _mirror_output_to_results(summary_path, os.path.join(compat_results_root, "comparison", "summary.json"))

        if manifest.get('methods'):
            method_statuses = [value.get("status") for value in manifest.get("methods", {}).values() if isinstance(value, dict)]
            if any(status == STEP_STATUS_SUCCESS for status in method_statuses):
                manifest['status'] = STEP_STATUS_SUCCESS
            elif method_statuses:
                manifest['status'] = STEP_STATUS_FAILED
            else:
                manifest['status'] = STEP_STATUS_SKIPPED
            manifest_path = _resolve_output_path(ps, "analysis", "final", "analysis_manifest.json")
            try:
                if analysis_artifacts:
                    manifest["artifacts"] = analysis_artifacts
                _write_json(manifest_path, manifest)
                _write_step_manifest(ps, "analysis", _build_step_manifest("analysis", manifest['status'], outputs={"methods": manifest.get('methods', {}), "artifacts": analysis_artifacts}, inputs={"mode": analysis_mode, "methods": methods}, error="; ".join(manifest.get("errors", []))))
            except Exception:
                pass
        else:
            manifest['status'] = STEP_STATUS_SKIPPED
            _write_step_manifest(ps, "analysis", _build_step_manifest("analysis", STEP_STATUS_SKIPPED, outputs={"methods": {}, "artifacts": analysis_artifacts}, inputs={"mode": analysis_mode, "methods": methods}, error="; ".join(manifest.get("errors", []))))

        return True
    except Exception as exc:
        try:
            _write_step_manifest(project_state, "analysis", _build_step_manifest("analysis", STEP_STATUS_FAILED, inputs={"mode": analysis_mode}, outputs={}, error=str(exc)))
        except Exception:
            pass
        raise

def run_mafft(input_fasta, output_fasta, timeout_sec=300):
    if not os.path.exists(input_fasta) or os.path.getsize(input_fasta) == 0:
        return False, "Input FASTA for MAFFT not found or empty."
    cmd = [MAFFT_PATH, "--auto", input_fasta]
    try:
        with open(output_fasta, "w") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, timeout=timeout_sec)
    except subprocess.TimeoutExpired:
        return False, f"MAFFT timed out after {timeout_sec}s."
    except Exception as e:
        return False, f"MAFFT error: {e}"
    if not os.path.exists(output_fasta) or os.path.getsize(output_fasta) == 0:
        return False, "MAFFT failed: output file not created."
    return True, result.stderr


def run_fasttree(input_fasta, output_tree, timeout_sec=600):
    from Bio import SeqIO
    try:
        seqs = list(SeqIO.parse(input_fasta, "fasta"))
    except Exception as e:
        return False, f"Error reading FASTA for FastTree: {e}"
    if len(seqs) < 2:
        return False, f"FastTree requires at least 2 sequences, found {len(seqs)}."
    cmd = [FASTTREE_PATH, "-gtr", "-nt", input_fasta]
    try:
        with open(output_tree, "w") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, timeout=timeout_sec)
    except subprocess.TimeoutExpired:
        return False, f"FastTree timed out after {timeout_sec}s."
    except Exception as e:
        return False, f"FastTree error: {e}"
    if not os.path.exists(output_tree) or os.path.getsize(output_tree) == 0:
        return False, "FastTree failed: output tree not created."
    return True, result.stderr


def run_raxml(input_fasta, output_dir, n_bootstrap=100, timeout_sec=3600):
    if not os.path.exists(input_fasta) or os.path.getsize(input_fasta) == 0:
        return False, "Input FASTA for RAxML-NG not found or empty."
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, "raxml")
    cmd = [RAXML_PATH, "--all", "--msa", input_fasta, "--model", "GTR+G",
           "--prefix", prefix, "--bs-trees", str(n_bootstrap), "--redo"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
    except subprocess.TimeoutExpired:
        return False, f"RAxML-NG timed out after {timeout_sec}s."
    except Exception as e:
        return False, f"RAxML-NG error: {e}"
    best_tree = prefix + ".raxml.bestTree"
    if not os.path.exists(best_tree) or os.path.getsize(best_tree) == 0:
        return False, result.stdout + "\n" + result.stderr
    return True, result.stdout + "\n" + result.stderr

def uniquify_fasta(in_fasta, out_fasta):
    from Bio import SeqIO
    seen = {}
    with open(out_fasta, "w") as out_handle:
        for i, rec in enumerate(SeqIO.parse(in_fasta, "fasta")):
            base_id = rec.id
            if base_id in seen:
                seen[base_id] += 1
                rec.id = f"{base_id}_{seen[base_id]}"
                rec.description = ""
            else:
                seen[base_id] = 1
            SeqIO.write(rec, out_handle, "fasta")

def filter_empty_aligned_fasta(in_fasta, out_fasta, min_real=10):
    from Bio import SeqIO
    records = []
    for rec in SeqIO.parse(in_fasta, "fasta"):
        real_bases = sum(1 for b in str(rec.seq) if b not in "-.")
        if real_bases >= min_real:
            records.append(rec)
    SeqIO.write(records, out_fasta, "fasta")

def sanitize_ids(in_fasta, out_fasta):
    from Bio import SeqIO
    import re
    with open(out_fasta, "w") as out_handle:
        for rec in SeqIO.parse(in_fasta, "fasta"):
            safe_id = re.sub(r'[^A-Za-z0-9_.-]', '_', rec.id)
            rec.id = safe_id
            rec.name = safe_id
            rec.description = ""
            SeqIO.write(rec, out_handle, "fasta")


def run_trimmomatic(input_file, output_file, trimmomatic_jar, adapters_file):
    cmd = [
        "java", "-jar", trimmomatic_jar, "SE", "-phred33",
        input_file, output_file,
        f"ILLUMINACLIP:{adapters_file}:2:30:10",
        "SLIDINGWINDOW:4:20", "MINLEN:50"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        result.stderr += "\n[ERROR] The trimmed.fastq file was not created or is empty. Check the input and Trimmomatic parameters."
    return result.stdout, result.stderr

def run_trimmomatic_pe(input_r1, input_r2, output_r1_paired, output_r1_unpaired, output_r2_paired, output_r2_unpaired, trimmomatic_jar, adapters_file):
    cmd = [
        "java", "-jar", trimmomatic_jar, "PE", "-phred33",
        input_r1, input_r2,
        output_r1_paired, output_r1_unpaired,
        output_r2_paired, output_r2_unpaired,
        f"ILLUMINACLIP:{adapters_file}:2:30:10",
        "SLIDINGWINDOW:4:20", "MINLEN:50"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout, result.stderr

def run_fastqc(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    cmd = [FASTQC_PATH, input_file, "-o", output_dir]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout, result.stderr

def run_pear(forward_file, reverse_file, output_dir, pear_path):
    cmd = [
        pear_path,
        "-f", forward_file,
        "-r", reverse_file,
        "-o", os.path.join(output_dir, "merged")
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout, result.stderr

def biopython_to_cytoscape_elements(tree):
    elements = []
    node_counter = 0

    def _clade_to_elements_recursive(clade, parent_id=None):
        nonlocal node_counter, elements
        node_counter += 1
        
        clade_name_safe = str(clade.name).replace('.', '_').replace(':', '_') if clade.name else f'internal_{node_counter}'
        current_id = f"node_{node_counter}_{clade_name_safe}"

        label = ""
        if clade.name:
            label = clade.name
        elif hasattr(clade, 'confidence') and clade.confidence is not None:
            try:
                label = str(round(float(clade.confidence), 2))
            except (ValueError, TypeError):
                label = str(clade.confidence) 
        
        elements.append({'data': {'id': current_id, 'label': label}})

        if parent_id:
            elements.append({'data': {'source': parent_id, 'target': current_id}})

        if clade.clades: 
            for child_clade in clade.clades:
                _clade_to_elements_recursive(child_clade, current_id)

    if tree and tree.root:
        _clade_to_elements_recursive(tree.root)
    return elements

@app.callback(
    [Output("project-status", "children"), Output("project-state", "data", allow_duplicate=True), Output("samples-store", "data", allow_duplicate=True)],
    Input("create-project-btn", "n_clicks"),
    State("project-name", "value"),
    State("project-state", "data"),
    prevent_initial_call=True
)
def create_or_select_project(n_clicks, project_name, project_state):
    if not n_clicks:
        return dash.no_update, dash.no_update, dash.no_update
    if not project_name or not str(project_name).strip():
        return html.Div("Please provide a project name before continuing.", style={"color": "orange"}), dash.no_update, dash.no_update

    pdir = _ensure_project(project_name)
    state = (project_state or {}).copy()
    state.update({"project_name": _safe_name(project_name), "project_dir": pdir})
    return html.Div(f"Active project set: {os.path.basename(pdir)}. Now upload ZIP/RAR dataset.", style={"color": "green"}), state, dash.no_update


@app.callback(
    Output("existing-projects", "options"),
    [Input("theme-root", "className"), Input("project-status", "children"), Input("zip-upload-status", "children")],
)
def refresh_existing_projects(_, __, ___):
    return _list_existing_projects()


@app.callback(
    [
        Output("project-load-status", "children"),
        Output("project-state", "data", allow_duplicate=True),
        Output("samples-store", "data", allow_duplicate=True),
        Output("active-sample", "data", allow_duplicate=True),
        Output("run-config", "data", allow_duplicate=True),
        Output("uploaded-files", "data", allow_duplicate=True),
        Output("run-mode", "value", allow_duplicate=True),
        Output("run-steps", "value", allow_duplicate=True),
        Output("project-name", "value", allow_duplicate=True),
        Output("workflow-scope-samples", "data", allow_duplicate=True),
    ],
    Input("load-project-btn", "n_clicks"),
    State("existing-projects", "value"),
    prevent_initial_call=True
)
def load_existing_project(n_clicks, selected_project):
    if not n_clicks or not selected_project:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    project_dir = _project_dir(selected_project)
    if not os.path.exists(project_dir):
        return html.Div("Selected project not found on disk.", style={"color": "red"}), dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    session_state, run_config, uploaded_files, active_sample, session_samples, session_scope = _load_project_session(project_dir)
    latest_run = _latest_run_dir(project_dir)
    dataset_dir = os.path.join(project_dir, "uploads", "dataset")
    indexed_samples = _index_samples_from_dir(dataset_dir) if os.path.exists(dataset_dir) else {}
    samples = indexed_samples or session_samples

    restored_state = {
        "project_name": _safe_name(selected_project),
        "project_dir": project_dir,
        "dataset_dir": dataset_dir if os.path.exists(dataset_dir) else "",
        "run_dir": (session_state or {}).get("run_dir") or latest_run or ""
    }

    mode = (run_config or {}).get("mode") or "complete"
    steps = (run_config or {}).get("steps") or WORKFLOW_STEP_ORDER.copy()
    if isinstance(active_sample, list):
        active_sample = [sid for sid in active_sample if sid in samples]
        if not active_sample:
            active_sample = "__ALL__"
    elif active_sample not in ("__ALL__",) and active_sample not in samples:
        active_sample = "__ALL__"

    # Ensure workflow_scope has the correct structure
    workflow_scope = session_scope or {}
    if not workflow_scope.get("samples_selected"):
        workflow_scope["samples_selected"] = []
    if "pe_pairs" not in workflow_scope:
        workflow_scope["pe_pairs"] = {}

    msg = html.Div(f"Project loaded: {selected_project}. Session restored.", style={"color": "green"})
    return msg, restored_state, samples, active_sample or "__ALL__", run_config or {}, uploaded_files or {"preprocess_mode": "se"}, mode, steps, selected_project, workflow_scope

@app.callback(
    [Output("zip-upload-status", "children"), Output("samples-store", "data", allow_duplicate=True), Output("project-state", "data", allow_duplicate=True), Output("active-sample", "data", allow_duplicate=True)],
    Input("project-zip-upload", "contents"),
    State("project-zip-upload", "filename"),
    State("project-state", "data"),
    prevent_initial_call=True
)
def handle_zip_upload(contents, filename, project_state):
    try:
        if not contents or not filename:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update
        state = (project_state or {}).copy()
        pdir = state.get("project_dir")
        if not pdir:
            return html.Div("Select/Create a project first, then upload the dataset archive.", style={"color": "orange"}), dash.no_update, dash.no_update, dash.no_update

        project_name = state.get("project_name") or os.path.basename(pdir)
        uploads_dir = os.path.join(pdir, "uploads")
        os.makedirs(uploads_dir, exist_ok=True)
        safe_archive_name = _safe_name(filename)
        archive_path = os.path.join(uploads_dir, safe_archive_name)

        try:
            _, content_string = contents.split(",", 1)
            decoded = base64.b64decode(content_string)
        except Exception as e:
            return html.Div(f"Invalid upload payload: {e}", style={"color": "red"}), dash.no_update, dash.no_update, dash.no_update

        with open(archive_path, "wb") as f:
            f.write(decoded)

        dataset_dir = os.path.join(uploads_dir, "dataset")
        if os.path.exists(dataset_dir):
            shutil.rmtree(dataset_dir)
        try:
            _safe_extract_archive(archive_path, dataset_dir)
        except Exception as e:
            return html.Div(f"Archive extraction error: {e}", style={"color": "red"}), dash.no_update, dash.no_update, dash.no_update

        samples = _index_samples_from_dir(dataset_dir)
        active_sample = "__ALL__"
        state.update({"project_name": _safe_name(project_name), "project_dir": pdir, "dataset_dir": dataset_dir})
        return html.Div(f"Dataset imported into project '{_safe_name(project_name)}': {filename} ({len(samples)} samples indexed)", style={"color": "green"}), samples, state, active_sample
    except Exception as e:
        return html.Div(f"An unexpected error occurred during ZIP upload: {e}", style={"color": "red"}), dash.no_update, dash.no_update, dash.no_update

@app.callback(
    Output("sample-selector", "options"),
    Input("samples-store", "data")
)
def update_sample_selector(samples_store):
    samples = samples_store or {}
    options = [{"label": "All samples (comparison mode)", "value": "__ALL__"}]
    options.extend(_group_samples_for_options(samples))
    return options


@app.callback(
    Output('samples-store', 'data', allow_duplicate=True),
    Input('samples-store', 'data'),
    prevent_initial_call=True
)
def normalize_samples_store_callback(samples_store):
    # Normalize and merge paired samples; avoid infinite update loop by comparing
    normalized = _normalize_samples_store(samples_store or {})
    # If identical, do not update
    if normalized == (samples_store or {}):
        return dash.no_update
    return normalized

@app.callback(
    Output("active-sample", "data", allow_duplicate=True),
    Input("sample-selector", "value"),
    prevent_initial_call=True
)
def set_active_sample(sample_id):
    return _normalize_active_sample_value(sample_id)


@app.callback(
    Output("sample-selector", "value"),
    Input("active-sample", "data"),
    prevent_initial_call=False
)
def sync_sample_selector_value(active_sample):
    if active_sample == "__ALL__" or not active_sample:
        return ["__ALL__"]
    if isinstance(active_sample, list):
        return active_sample
    return [active_sample]


@app.callback(
    Output("pipeline-state", "data"),
    [
        Input("project-state", "data"),
        Input("run-config", "data"),
        Input("uploaded-files", "data"),
        Input("active-sample", "data"),
        Input("samples-store", "data"),
        Input("workflow-scope-samples", "data"),
    ],
    prevent_initial_call=True
)
def persist_session_state(project_state, run_config, uploaded_files, active_sample, samples_store, workflow_scope):
    saved = _save_project_session(project_state, run_config, uploaded_files, active_sample, samples_store, workflow_scope)
    return {"saved": bool(saved), "saved_at": datetime.now().isoformat(timespec="seconds")}


@app.callback(
    Output("preprocess-sample-selector", "options"),
    Input("samples-store", "data"),
    Input("uploaded-files", "data"),
)
def update_preprocess_sample_options(samples_store, uploaded_files):
    samples = samples_store or {}
    # default top option
    options = [{"label": "All selected/available samples", "value": "__ALL__"}]
    grouped = _group_samples_for_options(samples_store)

    # Extract preprocess_mode from uploaded_files store (source of truth)
    preprocess_mode = (uploaded_files or {}).get('preprocess_mode', 'se')

    # filter grouped options depending on preprocess_mode
    if preprocess_mode == 'pe':
        # only include sample SIDs that have both r1 and r2
        for opt in grouped:
            sid = opt['value']
            info = samples.get(sid, {})
            paths = info.get('paths', {})
            if paths.get('r1') and paths.get('r2'):
                options.append(opt)
    else:
        # include all grouped options (SE entries will be labeled accordingly)
        options.extend(grouped)
    return options

@app.callback(
    [Output("run-status", "children"), Output("project-state", "data", allow_duplicate=True), Output("run-config", "data"), Output("workflow-scope-samples", "data", allow_duplicate=True)],
    Input("create-run-btn", "n_clicks"),
    State("project-state", "data"),
    State("run-mode", "value"),
    State("run-steps", "value"),
    State("project-name", "value"),
    State("workflow-scope-samples", "data"),
    State("samples-store", "data"),
    prevent_initial_call=True
)
def create_run(n_clicks, project_state, run_mode, run_steps, project_name, workflow_scope, samples_store):
    if not n_clicks:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update
    state = (project_state or {}).copy()
    pdir = state.get("project_dir") or _ensure_project(project_name or state.get("project_name") or "Project")

    steps = run_steps or []
    if run_mode == "complete":
        steps = ["qc", "preprocess", "analysis", "db", "phylo"]
    elif run_mode == "special" and len(steps) < 2:
        return html.Div("Special mode requires at least 2 steps.", style={"color": "orange"}), dash.no_update, dash.no_update, dash.no_update
    elif run_mode == "one" and len(steps) != 1:
        return html.Div("One mode requires exactly 1 step.", style={"color": "orange"}), dash.no_update, dash.no_update, dash.no_update

    rdir = _create_run(pdir)
    config = {
        "mode": run_mode,
        "steps": steps,
        "created_at": datetime.now().isoformat(timespec="seconds")
    }
    with open(os.path.join(rdir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    state.update({"project_name": state.get("project_name") or _safe_name(project_name), "project_dir": pdir, "run_dir": rdir})
    
    # Use the existing workflow_scope from the store; if empty, initialize
    if not workflow_scope or not workflow_scope.get("samples_selected"):
        workflow_scope = {
            "samples_selected": [],
            "selected_at": datetime.now().isoformat(timespec="seconds"),
            "pe_pairs": {}
        }

    _write_run_session_metadata(state, config, workflow_scope, samples_store, outputs={})
    
    return html.Div(f"Run created: {os.path.basename(rdir)}"), state, config, workflow_scope

# ----------------------------------------------------------------------
# NEW: Workflow Sample Scope Card and callbacks
# ----------------------------------------------------------------------

@app.callback(
    Output("workflow-scope-selector", "options"),
    Output("workflow-scope-selector", "value"),
    Output("workflow-scope-counter", "children"),
    Output("workflow-scope-samples", "data", allow_duplicate=True),
    Input("samples-store", "data"),
    Input("workflow-scope-search", "value"),
    Input("workflow-scope-select-all", "n_clicks"),
    Input("workflow-scope-clear", "n_clicks"),
    Input("workflow-scope-selector", "value"),
    State("workflow-scope-samples", "data"),
    prevent_initial_call=True  # <--- CORRIGIDO: adicionado prevent_initial_call=True
)
def update_workflow_scope(
    samples_store,
    search_text,
    select_all_clicks,
    clear_clicks,
    selected_values,
    current_scope
):
    ctx = dash.callback_context
    samples_store = samples_store or {}
    current_scope = current_scope or {}

    sample_ids = sorted(samples_store.keys())
    if search_text:
        q = search_text.lower()
        filtered = [
            sid for sid in sample_ids
            if q in sid.lower()
            or q in str((samples_store.get(sid, {}) or {}).get('display_name', '')).lower()
            or q in str((samples_store.get(sid, {}) or {}).get('display_label', '')).lower()
            or q in str((samples_store.get(sid, {}) or {}).get('internal_sample_id', '')).lower()
        ]
    else:
        filtered = sample_ids

    # Use grouped, human-friendly labels but keep option values as original SIDs
    options = _group_samples_for_options({k: samples_store.get(k, {}) for k in filtered})

    trigger = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else None
    if trigger == "workflow-scope-select-all":
        selected_values = filtered
    elif trigger == "workflow-scope-clear":
        selected_values = []

    selected_values = selected_values or []
    # Build pe_pairs mapping for quick lookup: sid -> member file list
    pe_pairs = {}
    for sid in selected_values:
        info = samples_store.get(sid, {})
        pe_pairs[sid] = {
            'members': info.get('members', []),
            'paths': info.get('paths', {})
        }
    scope_payload = {
        "samples_selected": selected_values,
        "selected_at": datetime.now().isoformat(timespec="seconds"),
        "pe_pairs": pe_pairs
    }
    counter = f"{len(selected_values)} samples selected"
    return options, selected_values, counter, scope_payload

# ----------------------------------------------------------------------

app.layout = html.Div(
    id="theme-root",
    className="theme-light",
    children=[
        dbc.Container([
            dcc.Store(id="project-state", data={}),
            dcc.Store(id="run-config", data={}),
            dcc.Store(id="samples-store", data={}),
            dcc.Store(id="pipeline-state", data={}),
            dcc.Store(id="uploaded-files", data={"preprocess_mode": "se"}),
            dcc.Store(id="active-sample", data="__ALL__"),
            dcc.Store(id="workflow-progress", data={
                "qc_done": False,
                "preprocess_done": False,
                "analysis_done": False,
                "db_done": False,
                "phylo_done": False,
                "results_done": False
            }),
            dcc.Store(id="temp-r1-path-store", data=None),
            dcc.Store(id="temp-r2-path-store", data=None),
            dcc.Store(id="temp-se-status-store", data=None),
            dcc.Store(id="workflow-scope-samples", data={
                "samples_selected": [],
                "selected_at": None,
                "pe_pairs": {}
            }),
            dcc.Store(id="run-completed-signal", data=False),
            html.Div([
                html.Div([
                    html.Img(src="/uploads/logo.png", style={"height": "70px", "marginRight": "8px", "verticalAlign": "middle", "display": "inline-block"}),
                    html.H1("PARAFISH", style={"display": "inline-block", "verticalAlign": "middle", "margin": "0"}),
                ], style={"display": "flex", "alignItems": "center", "gap": "0.5rem"}),
                dbc.Switch(
                    id="theme-toggle",
                    label="Dark mode",
                    value=False,
                    className="mb-0",
                    persistence=True,
                    persistence_type="local"
                )
            ], style={"display": "flex", "justifyContent": "space-between", "alignItems": "center", "marginTop": "10px", "marginBottom": "10px"}),

            dbc.Card([
                dbc.CardHeader("1) Setup: Project & Dataset"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            html.Label("Project name"),
                            dbc.Input(id="project-name", placeholder="Project name", type="text"),
                            dbc.Button("Create/Select Project", id="create-project-btn", color="primary", className="mt-2"),
                            html.Div(id="project-status", className="mt-2"),
                            html.Hr(),
                            html.Label("Load existing project"),
                            dcc.Dropdown(id="existing-projects", placeholder="Select an existing project", clearable=True),
                            dbc.Button("Load Project", id="load-project-btn", color="secondary", className="mt-2"),
                            html.Div(id="project-load-status", className="mt-2")
                        ], width=4),
                        dbc.Col([
                            html.Label("Upload dataset archive (ZIP/RAR)"),
                            dcc.Upload(
                                id="project-zip-upload",
                                children=html.Div(["Drag and Drop or ", html.A("Select ZIP/RAR")]),
                                style={"width": "100%", "height": "60px", "lineHeight": "60px", "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "5px", "textAlign": "center"},
                                multiple=False
                            ),
                            html.Div(id="zip-upload-status", className="mt-2")
                        ], width=4),
                        dbc.Col([
                            html.Label("Active sample (preview only)"),
                            dcc.Dropdown(id="sample-selector", placeholder="Select indexed sample(s)", clearable=False, value=["__ALL__"], multi=True),
                            html.Div(id="setup-gate-status", className="mt-2")
                        ], width=4),
                    ])
                ])
            ], className="mb-3"),

            dbc.Card([
                dbc.CardHeader("2) Workflow Mode & Steps"),
                dbc.CardBody([
                    html.Div(id="setup-required-message", className="mb-2"),
                    dbc.Row([
                        dbc.Col([
                            html.Label("Execution mode"),
                            dbc.RadioItems(
                                id="run-mode",
                                options=[
                                    {"label": "Complete", "value": "complete"},
                                    {"label": "Special", "value": "special"},
                                    {"label": "One", "value": "one"}
                                ],
                                value="complete",
                                inline=True
                            )
                        ], width=4),
                        dbc.Col([
                            html.Label("Select steps"),
                            dbc.Checklist(
                                id="run-steps",
                                options=[
                                    {"label": "QC", "value": "qc"},
                                    {"label": "Preprocess", "value": "preprocess"},
                                    {"label": "Analysis", "value": "analysis"},
                                    {"label": "Databases", "value": "db"},
                                    {"label": "Phylo", "value": "phylo"}
                                ],
                                value=["qc", "preprocess", "analysis", "db", "phylo"],
                                inline=True
                            ),
                            html.Div(id="run-validation", className="mt-2")
                        ], width=8),
                    ]),
                    html.Hr(),
                    dbc.Button("Create Run", id="create-run-btn", color="secondary"),
                    html.Div(id="run-status", className="mt-2")
                ])
            ], id="workflow-config-card", className="mb-3"),

            # NEW: Workflow Sample Scope Card
            dbc.Card([
                dbc.CardHeader(html.H5("3) Workflow Sample Scope")),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dbc.Input(
                                id="workflow-scope-search",
                                placeholder="Search samples...",
                                type="text"
                            )
                        ], md=6),
                        dbc.Col([
                            dbc.ButtonGroup([
                                dbc.Button("Select All", id="workflow-scope-select-all", color="primary", outline=True),
                                dbc.Button("Clear", id="workflow-scope-clear", color="secondary", outline=True)
                            ])
                        ], md=6, className="text-end")
                    ], className="mb-3"),
                    dcc.Checklist(
                        id="workflow-scope-selector",
                        options=[],
                        value=[],
                        style={
                            "maxHeight": "300px",
                            "overflowY": "auto",
                            "padding": "10px",
                            "border": "1px solid #ddd",
                            "borderRadius": "8px"
                        }
                    ),
                    html.Div(id="workflow-scope-counter", className="mt-3 text-muted")
                ])
            ], className="mb-3"),

            dbc.Card([
                dbc.CardHeader("4) Selected Steps (Accordion Workflow)"),
                dbc.CardBody([
                    html.Div(id="selected-steps-summary", className="mb-2"),
                    html.Div(
                        id="selected-steps-panel",
                        children=html.Div([
                            html.Div(id='upload-status-qc'),
                            html.Div(id='upload-status-analysis'),
                            html.Div(id='upload-status-databases'),
                            html.Div(id='upload-status-phylo'),
                            html.Div(id='upload-status-preprocess'),
                            html.Div(id='preprocess-uploader-container'),
                            dcc.Dropdown(id='preprocess-sample-selector'),
                            html.Div(id='qc-results'),
                            html.Div(id='qc-summary-results'),
                            html.Div(id='preprocess-results'),
                            html.Div(id='analysis-results'),
                            html.Div(id='combined-fasta-status'),
                            html.Div(id='db-info'),
                            html.Div(id='custom-db-status'),
                            dcc.Dropdown(id='custom-db-dropdown'),
                            html.Div(id='db-classify-results'),
                            html.Div(id='boldigger-results'),
                            html.Div(id='phylo-results')
                        ], style={"display": "none"})
                    ),
                    html.Hr(),
                    dbc.Button("Run Selected Workflow", id="run-pipeline-btn", color="primary"),
                    html.Div(id="export-status", className="mt-2")
                ])
            ], id="workflow-steps-card", className="mb-3"),

            dbc.Card([
                dbc.CardHeader("Results"),
                dbc.CardBody([
                    dcc.Loading(
                        id="pipeline-loading",
                        type="default",
                        children=html.Div(id="pipeline-results")
                    ),
                    html.Hr(),
                    html.Div(id="final-results")
                ])
            ])
        ], fluid=True)
    ]
)

# --- Multi-sample constants and helpers ---

def _next_sample_id(existing_ids):
    i = 1
    while True:
        sid = f"S{str(i).zfill(3)}"
        if sid not in existing_ids:
            return sid
        i += 1

def assign_samples_from_files(file_list, current_samples=None, mode_hint=None):
    samples = (current_samples or {}).copy()
    existing_ids = set(samples.keys())

    def _display_label(internal_sample_id, marker, mode):
        label = str(internal_sample_id or "sample")
        if marker:
            label = f"{label} ({marker})"
        return f"{label} ({mode})"

    for fp in file_list or []:
        if not fp:
            continue
        parsed = parse_fastq_name(fp)
        base = os.path.basename(fp)
        pair_key = parsed.internal_sample_id
        mode = str(mode_hint or parsed.mode or "SE").upper()

        # try to find existing sample by exact pair_key match
        sid = None
        for k, v in samples.items():
            if v.get('pair_key') == pair_key:
                sid = k
                break

        if not sid:
            sid = _next_sample_id(existing_ids)
            existing_ids.add(sid)
            samples[sid] = {
                'display_name': _display_label(pair_key, parsed.marker, mode),
                'display_label': _display_label(pair_key, parsed.marker, mode),
                'stem': os.path.splitext(base)[0],
                'pair_key': pair_key,
                'internal_sample_id': pair_key,
                'kind': mode,
                'mode': mode,
                'marker': parsed.marker,
                'members': [parsed.display_label],
                'paths': {},
                'source_files': [fp],
                'outputs': {},
                'meta': {}
            }
        else:
            # update display name to include both members if pairing
            if parsed.display_label not in samples[sid].get('members', []):
                samples[sid].setdefault('members', []).append(parsed.display_label)
                samples[sid]['display_name'] = _display_label(pair_key, samples[sid].get('marker') or parsed.marker, samples[sid].get('mode', mode))
                samples[sid]['display_label'] = samples[sid]['display_name']
            samples[sid].setdefault('source_files', [])
            if fp not in samples[sid]['source_files']:
                samples[sid]['source_files'].append(fp)
            if not samples[sid].get('internal_sample_id'):
                samples[sid]['internal_sample_id'] = pair_key
            if not samples[sid].get('marker') and parsed.marker:
                samples[sid]['marker'] = parsed.marker

        # assign file paths
        if parsed.mode == 'PE' and parsed.mate == 1:
            samples[sid]['paths']['r1'] = fp
            samples[sid]['kind'] = 'PE'
            samples[sid]['mode'] = 'PE'
        elif parsed.mode == 'PE' and parsed.mate == 2:
            samples[sid]['paths']['r2'] = fp
            samples[sid]['kind'] = 'PE'
            samples[sid]['mode'] = 'PE'
        else:
            samples[sid]['paths']['se'] = fp
            # if a SE file exists alongside r1/r2, prefer marking as SE unless r1/r2 present
            if not samples[sid]['paths'].get('r1') and not samples[sid]['paths'].get('r2'):
                samples[sid]['kind'] = 'SE'
                samples[sid]['mode'] = 'SE'

        # ensure display_name present
        if 'display_name' not in samples[sid] or not samples[sid]['display_name']:
            marker = samples[sid].get('marker')
            samples[sid]['display_name'] = _display_label(samples[sid].get('internal_sample_id') or pair_key, marker, samples[sid].get('mode', 'SE'))
            samples[sid]['display_label'] = samples[sid]['display_name']

    return samples

def write_combined_fasta(sample_fastas, out_path):
    from Bio import SeqIO
    with open(out_path, 'w') as out_f:
        for f in sample_fastas:
            if not f or not os.path.exists(f) or os.path.getsize(f) == 0:
                continue
            for rec in SeqIO.parse(f, 'fasta'):
                SeqIO.write(rec, out_f, 'fasta')
    return os.path.exists(out_path) and os.path.getsize(out_path) > 0


WORKFLOW_STEP_ORDER = ["qc", "preprocess", "analysis", "db", "phylo"]
WORKFLOW_STEP_LABELS = {
    "qc": "QC",
    "preprocess": "Preprocess",
    "analysis": "Analysis",
    "db": "Databases",
    "phylo": "Phylogenetics"
}


def _is_setup_ready(project_state, samples_store):
    state = project_state or {}
    has_project = bool(state.get("project_dir"))
    has_dataset = bool(state.get("dataset_dir")) or bool(samples_store)
    return has_project and has_dataset


def _normalize_steps(run_mode, run_steps):
    selected = []
    for step in (run_steps or []):
        if step in WORKFLOW_STEP_ORDER and step not in selected:
            selected.append(step)

    if run_mode == "complete":
        return WORKFLOW_STEP_ORDER.copy()
    if run_mode == "one":
        return selected[:1]
    return selected


def _analysis_mode_to_methods(mode_value):
    mode = str(mode_value or "vsearch_only").strip().lower()
    if mode == "dada2_only":
        return ["dada2"]
    if mode == "both":
        return ["vsearch", "dada2"]
    return ["vsearch"]


def _render_step_content(step_key):
    if step_key == "qc":
        return html.Div([
            du.Upload(
                id='du-upload-qc',
                text='Upload FASTQ for QC (optional)',
                max_files=1,
                max_file_size=100000,
                filetypes=['fastq.gz', 'gz', 'fastq'],
            ),
            html.Div(id='upload-status-qc', className="mt-2"),
            html.P("Execution is managed by 'Run Selected Workflow'.", style={"color": "#6c757d"}),
            dbc.Button("Run Quality Control", id="run-qc-btn", color="primary", style={"display": "none"}),
            dbc.Button("Build QC Summary", id="build-qc-summary-btn", color="secondary", className="ms-2", style={"display": "none"}),
            html.Div(id="qc-summary-results", className="mt-3"),
            dcc.Loading(type="default", children=html.Div(id="qc-results", className="mt-3")),
        ])

    if step_key == "preprocess":
        return html.Div([
            html.Label("Preprocess mode"),
            dbc.RadioItems(
                id="preprocess-mode",
                options=[
                    {"label": "Single-end", "value": "se"},
                    {"label": "Paired-end", "value": "pe"},
                ],
                value="se",
                inline=True,
                className="mb-2"
            ),
            dbc.Button("Confirm mode", id="preprocess-confirm-btn", color="info", className="mb-3"),
            html.Div([
                html.Label("Preprocess sample (for PE choose a paired sample entry, for SE choose a single-end sample)"),
                dcc.Dropdown(id="preprocess-sample-selector", clearable=True, value="", className="mb-2"),
            ], id="preprocess-pe-overrides-container"),
            html.Div(id="preprocess-uploader-container"),
            html.Div(id='upload-status-preprocess', className="mt-2"),
            html.P("Execution is managed by 'Run Selected Workflow'.", style={"color": "#6c757d"}),
            dbc.Button("Run Trimmomatic", id="run-trim-btn", color="primary", className="me-2", style={"display": "none"}),
            dbc.Button("Run PEAR", id="run-pear-btn", color="secondary", style={"display": "none"}),
            dcc.Loading(type="default", children=html.Div(id="preprocess-results", className="mt-3")),
        ])

    if step_key == "analysis":
        return html.Div([
            du.Upload(
                id='du-upload-analysis',
                text='Upload FASTA for Analysis (optional)',
                max_files=1,
                max_file_size=100000,
                filetypes=['fasta', 'fa'],
            ),
            html.Div(id='upload-status-analysis', className="mt-2"),
            dbc.RadioItems(
                id="analysis-method",
                options=[
                    {"label": "vsearch only (OTUs)", "value": "vsearch_only"},
                    {"label": "DADA2 only (ASVs)", "value": "dada2_only"},
                    {"label": "both (independent OTU + ASV)", "value": "both"},
                ],
                value="vsearch_only",
                inline=True
            ),
            html.Div([
                html.Label("Identity threshold (OTU clustering)"),
                dcc.Slider(
                    id="vsearch-identity",
                    min=0.8, max=1.0, step=0.01, value=0.97,
                    marks={0.8: "80%", 0.85: "85%", 0.9: "90%", 0.95: "95%", 0.97: "97%", 1.0: "100%"},
                    tooltip={"placement": "bottom", "always_visible": True}
                )
            ], style={"marginTop": "1em", "marginBottom": "1em"}),
            html.P("Execution is managed by 'Run Selected Workflow'.", style={"color": "#6c757d"}),
            dbc.Button("Run Analysis", id="run-analysis-btn", color="primary", style={"display": "none"}),
            html.Hr(),
            dbc.Button("Build Combined FASTA", id="build-combined-fasta-btn", color="secondary", style={"display": "none"}),
            html.Div(id='combined-fasta-status', className='mt-2'),
            dcc.Loading(type="default", children=html.Div(id="analysis-results", className="mt-3")),
        ])

    if step_key == "db":
        return html.Div([
            du.Upload(
                id='du-upload-databases',
                text='Upload FASTA for Databases step (optional)',
                max_files=1,
                max_file_size=100000,
                filetypes=['fasta', 'fa'],
            ),
            html.Div(id='upload-status-databases', className="mt-2"),
            dbc.Checklist(
                id="db-selection",
                options=[
                    {"label": "NCBI GenBank", "value": "ncbi"},
                    {"label": "SILVA", "value": "silva"},
                    {"label": "PR2", "value": "pr2"},
                    {"label": "BOLD (via BOLDigger)", "value": "bold"},
                    {"label": "Custom DB", "value": "custom"}
                ],
                value=[],
                inline=True
            ),
            html.Div(id="db-info", className="mt-2"),
            html.Div([
                html.Label("BLAST max hits"),
                dcc.Input(id="blast-max-hits", type="number", min=1, max=50, step=1, value=5, style={"width": "120px"}),
                html.Label("BLAST e-value"),
                dcc.Input(id="blast-evalue", type="text", value="1e-10", style={"width": "120px"}),
            ], style={"display": "flex", "alignItems": "center", "gap": "0.5rem", "marginTop": "0.75em"}),
            html.Hr(),
            html.H6("Custom DB"),
            dcc.Upload(
                id="custom-db-upload",
                children=html.Div(["Drag and Drop or ", html.A("Select FASTA/ZIP")]),
                style={"width": "100%", "height": "60px", "lineHeight": "60px", "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "5px", "textAlign": "center"},
                multiple=False
            ),
            html.Div(id="custom-db-status", className="mt-2"),
            html.Label("Available custom DBs"),
            dcc.Dropdown(id="custom-db-dropdown", placeholder="Select a custom DB", clearable=True),
            html.Div(
                id="boldigger-options-container",
                children=[
                    html.Hr(),
                    html.H6("BOLDigger3 Options"),
                    html.Div([
                        html.Label("BOLD Database (1-8)"),
                        dcc.Dropdown(
                            id="bold-db-dropdown",
                            options=[
                                {"label": f"{i}: {name}", "value": i} for i, name in [
                                    (1, "ANIMAL LIBRARY (PUBLIC)"),
                                    (2, "ANIMAL SPECIES-LEVEL LIBRARY (PUBLIC + PRIVATE)"),
                                    (3, "ANIMAL LIBRARY (PUBLIC+PRIVATE)"),
                                    (4, "VALIDATED CANADIAN ARTHROPOD LIBRARY"),
                                    (5, "PLANT LIBRARY (PUBLIC)"),
                                    (6, "FUNGI LIBRARY (PUBLIC)"),
                                    (7, "ANIMAL SECONDARY MARKERS (PUBLIC)"),
                                    (8, "VALIDATED ANIMAL RED LIST LIBRARY"),
                                ]
                            ],
                            value=1,
                            clearable=False
                        ),
                        html.Label("BOLD Operating Mode (1-3)"),
                        dcc.Dropdown(
                            id="bold-mode-dropdown",
                            options=[
                                {"label": "1: Rapid Species Search", "value": 1},
                                {"label": "2: Genus and Species Search", "value": 2},
                                {"label": "3: Exhaustive Search", "value": 3}
                            ],
                            value=1,
                            clearable=False
                        ),
                        html.Label("Custom Thresholds (Species, Genus, Family, Order, Class)"),
                        dcc.Input(id="bold-thresholds", type="text", placeholder="Ex: 99 97 90 85 80", style={"width": "100%"})
                    ], style={"display": "grid", "gap": "0.5rem"}),
                ]
            ),
            html.Hr(),
            html.P("Execution is managed by 'Run Selected Workflow'.", style={"color": "#6c757d"}),
            dbc.Button("Run BOLDigger3", id="run-boldigger-btn", color="primary", style={"display": "none"}),
            dcc.Loading(type="default", children=html.Div(id="boldigger-results", className="mt-3")),
            html.Hr(),
            dbc.Button("Run DB Classification", id="run-db-classify-btn", color="secondary", style={"display": "none"}),
            dcc.Loading(type="default", children=html.Div(id="db-classify-results", className="mt-3")),
        ])

    if step_key == "phylo":
        return html.Div([
            du.Upload(
                id='du-upload-phylo',
                text='Upload FASTA for Phylogenetics (optional)',
                max_files=1,
                max_file_size=100000,
                filetypes=['fasta', 'fa'],
            ),
            html.Div(id='upload-status-phylo', className="mt-2"),
            html.Label("Phylogeny input source"),
            dcc.Dropdown(
                id='phylo-input-choice',
                options=[
                    {"label": "Auto (prefer OTUs then ASVs)", "value": "auto"},
                    {"label": "Use OTUs (VSEARCH)", "value": "otu"},
                    {"label": "Use ASVs (DADA2)", "value": "asv"}
                ],
                value='auto',
                clearable=False,
                style={"width": "320px"}
            ),
            html.Label("Number of bootstraps for RAxML-NG"),
            dcc.Input(
                id="raxml-bootstrap-input",
                type="number",
                min=10,
                max=1000,
                step=10,
                value=100,
                style={"width": "140px"}
            ),
            html.P("Execution is managed by 'Run Selected Workflow'.", style={"color": "#6c757d"}),
            dbc.Button("Run Phylogenetic Analysis", id="run-phylo-btn", color="primary", className="ms-2", style={"display": "none"}),
            dcc.Loading(type="default", children=html.Div(id="phylo-results", className="mt-3")),
        ])

    return html.Div()


def _run_preprocess_auto(uploaded_files, project_state, samples_store, workflow_scope, explicit_pe_sample=None):
    data = uploaded_files or {}
    mode = data.get("preprocess_mode", "se")
    target_samples = _iter_workflow_samples(samples_store, workflow_scope)
    if not target_samples:
        return html.Div("No samples available for preprocessing.", style={"color": "orange"})

    results_blocks = []
    base_dir = _step_base(project_state, "preprocess")
    
    # Use explicit sample if provided, otherwise fall back to store paths
    final_pe_sample = explicit_pe_sample
    explicit_pe_consumed = False

    for sid, sample in target_samples:
        paths = (sample or {}).get("paths", {})

        if mode == "pe":
            r1, r2 = None, None
            if final_pe_sample:
                if not explicit_pe_consumed:
                    sample_info = (samples_store or {}).get(final_pe_sample, {})
                    r1, r2 = sample_info.get('paths', {}).get('r1'), sample_info.get('paths', {}).get('r2')
                    explicit_pe_consumed = True
                else:
                    continue
            else:
                # fall back to stored paths per sample
                r1, r2 = paths.get("r1"), paths.get("r2")

            if not r1 or not r2:
                results_blocks.append(html.Div(f"{sid}: missing R1/R2 for paired-end preprocess.", style={"color": "orange"}))
                continue

            base1 = os.path.splitext(os.path.basename(r1))[0]
            base2 = os.path.splitext(os.path.basename(r2))[0]
            out_r1_paired = os.path.join(base_dir, f"{base1}_trimmed_R1_paired.fastq")
            out_r2_paired = os.path.join(base_dir, f"{base2}_trimmed_R2_paired.fastq")
            out_r1_unpaired = os.path.join(base_dir, f"{base1}_trimmed_R1_unpaired.fastq")
            out_r2_unpaired = os.path.join(base_dir, f"{base2}_trimmed_R2_unpaired.fastq")

            trim_stdout, trim_stderr = run_trimmomatic_pe(
                r1, r2, out_r1_paired, out_r1_unpaired, out_r2_paired, out_r2_unpaired,
                TRIMMOMATIC_JAR, TRIMMOMATIC_ADAPTERS_PE
            )
            pear_dir = _resolve_output_path(project_state, "preprocess", "intermediate", f"pear_results/{base1}_{base2}")
            os.makedirs(pear_dir, exist_ok=True)
            pear_stdout, pear_stderr = run_pear(out_r1_paired, out_r2_paired, pear_dir, PEAR_PATH)
            merged = os.path.join(pear_dir, "merged.assembled.fastq")

            results_blocks.append(html.Div([
                html.H6(f"{sid} (PE)"),
                html.A("Download merged PEAR output", href=_download_href(merged), target="_blank", download=os.path.basename(merged)),
                dbc.Accordion([
                    dbc.AccordionItem([html.H6("Trimmomatic stdout"), html.Pre(trim_stdout), html.H6("Trimmomatic stderr"), html.Pre(trim_stderr)], title="Log Trimmomatic"),
                    dbc.AccordionItem([html.H6("PEAR stdout"), html.Pre(pear_stdout), html.H6("PEAR stderr"), html.Pre(pear_stderr)], title="Log PEAR"),
                ], start_collapsed=True, always_open=False, style={"marginTop": "0.5em", "marginBottom": "0.5em"})
            ]))
            continue

        se_file = paths.get("se")
        if not se_file:
            se_list = data.get("preprocess", [])
            if se_list and isinstance(se_list, list):
                se_file = se_list[0]
        if not se_file:
            results_blocks.append(html.Div(f"{sid}: missing single-end FASTQ.", style={"color": "orange"}))
            continue

        base = os.path.splitext(os.path.basename(se_file))[0]
        output_file = os.path.join(base_dir, f"{base}_trimmed.fastq")
        trim_stdout, trim_stderr = run_trimmomatic(se_file, output_file, TRIMMOMATIC_JAR, TRIMMOMATIC_ADAPTERS)
        results_blocks.append(html.Div([
            html.H6(f"{sid} (SE)"),
            html.A("Download trimmed output", href=_download_href(output_file), target="_blank", download=os.path.basename(output_file)),
            dbc.Accordion([
                dbc.AccordionItem([html.H6("Trimmomatic stdout"), html.Pre(trim_stdout), html.H6("Trimmomatic stderr"), html.Pre(trim_stderr)], title="Log Trimmomatic"),
            ], start_collapsed=True, always_open=False, style={"marginTop": "0.5em", "marginBottom": "0.5em"})
        ]))

    if not results_blocks:
        return html.Div("No preprocess outputs were generated.", style={"color": "orange"})
    return html.Div(results_blocks)


@app.callback(
    [Output("setup-gate-status", "children"), Output("setup-required-message", "children")],
    [Input("project-state", "data"), Input("samples-store", "data")],
)
def setup_gate_status(project_state, samples_store):
    if _is_setup_ready(project_state, samples_store):
        ok = html.Span("Setup complete. You can now choose mode and steps.", style={"color": "green"})
        return ok, ok
    msg = html.Span("Complete setup first: create/select project and upload dataset archive (ZIP/RAR).", style={"color": "orange"})
    return msg, msg


@app.callback(
    [Output("project-zip-upload", "disabled"), Output("project-zip-upload", "style")],
    Input("project-state", "data")
)
def toggle_zip_upload_enabled(project_state):
    base_style = {
        "width": "100%",
        "height": "60px",
        "lineHeight": "60px",
        "borderWidth": "1px",
        "borderStyle": "dashed",
        "borderRadius": "5px",
        "textAlign": "center"
    }
    if (project_state or {}).get("project_dir"):
        return False, base_style
    disabled_style = base_style.copy()
    disabled_style.update({"opacity": 0.6, "cursor": "not-allowed"})
    return True, disabled_style


@app.callback(
    [Output("workflow-config-card", "style"), Output("workflow-steps-card", "style")],
    [Input("project-state", "data"), Input("samples-store", "data")],
)
def toggle_workflow_cards(project_state, samples_store):
    if _is_setup_ready(project_state, samples_store):
        return {}, {}
    return {"display": "none"}, {"display": "none"}


@app.callback(
    [Output("run-steps", "value"), Output("selected-steps-summary", "children"), Output("run-validation", "children")],
    [Input("run-mode", "value"), Input("run-steps", "value")],
    [State("project-state", "data"), State("samples-store", "data")],
    prevent_initial_call=True,
)
def validate_and_summarize_steps(run_mode, run_steps, project_state, samples_store):
    normalized = _normalize_steps(run_mode, run_steps)

    if not _is_setup_ready(project_state, samples_store):
        return normalized, html.Div("Workflow locked until setup is complete."), html.Div("Waiting for setup...", style={"color": "orange"})

    validation = html.Div("Selection is valid.", style={"color": "green"})
    if run_mode == "one" and len(normalized) != 1:
        validation = html.Div("Mode ONE requires exactly 1 step.", style={"color": "orange"})
    elif run_mode == "special" and len(normalized) < 2:
        validation = html.Div("Mode SPECIAL requires at least 2 steps.", style={"color": "orange"})
    elif run_mode == "complete":
        validation = html.Div("Mode COMPLETE includes all steps.", style={"color": "green"})

    if normalized:
        summary = html.Div([
            html.Span("Selected steps: "),
            html.Strong(" → ".join([WORKFLOW_STEP_LABELS[s] for s in normalized]))
        ])
    else:
        summary = html.Div("No steps selected yet.")

    return normalized, summary, validation


@app.callback(
    Output("selected-steps-panel", "children"),
    [Input("run-mode", "value"), Input("run-steps", "value"), Input("project-state", "data"), Input("samples-store", "data")],
)
def render_selected_steps_panel(run_mode, run_steps, project_state, samples_store):
    setup_ready = _is_setup_ready(project_state, samples_store)
    selected = _normalize_steps(run_mode, run_steps)

    items = []
    for idx, step in enumerate(WORKFLOW_STEP_ORDER, start=1):
        is_selected = step in selected

        if not setup_ready:
            body = html.Div([
                html.Div("Locked: complete setup (project + dataset) to enable this stage.", style={"color": "orange", "marginBottom": "0.5em"}),
                html.Div(_render_step_content(step), style={"display": "none"})
            ])
        elif not is_selected:
            body = html.Div([
                html.Div("Not selected in current mode.", style={"color": "#6c757d", "marginBottom": "0.5em"}),
                html.Div(_render_step_content(step), style={"display": "none"})
            ])
        else:
            body = _render_step_content(step)

        title_suffix = ""
        if not setup_ready:
            title_suffix = " (locked)"
        elif not is_selected:
            title_suffix = " (not selected)"

        items.append(
            dbc.AccordionItem(
                body,
                title=f"{idx}) {WORKFLOW_STEP_LABELS[step]}{title_suffix}"
            )
        )

    return dbc.Accordion(items, start_collapsed=True, always_open=True)


@app.callback(
    [Output("pipeline-results", "children"), Output("export-status", "children"), Output("run-completed-signal", "data")],
    Input("run-pipeline-btn", "n_clicks"),
    [State("run-config", "data"), State("run-mode", "value"), State("run-steps", "value"),
     State("uploaded-files", "data"), State("project-state", "data"), State("samples-store", "data"),
     State("workflow-scope-samples", "data"),
     State("preprocess-sample-selector", "value"),
     State("db-selection", "value"), State("custom-db-dropdown", "value"),
     State("blast-max-hits", "value"), State("blast-evalue", "value"),
     State("bold-db-dropdown", "value"), State("bold-mode-dropdown", "value"),
     State("bold-thresholds", "value"),
     State("analysis-method", "value"),
     State("vsearch-identity", "value"),
     State("raxml-bootstrap-input", "value"),
     State("phylo-input-choice", "value")],
    prevent_initial_call=True
)
def run_selected_workflow(
     n_clicks, run_config, run_mode, run_steps,
     uploaded_files, project_state, samples_store, workflow_scope,
     preprocess_sample,
     db_selection, custom_db, blast_max_hits, blast_evalue,
     bold_db, bold_mode, bold_thresholds,
     analysis_method, vsearch_identity, raxml_bootstrap, phylo_input_choice
 ):
    try:
        if not n_clicks:
            return dash.no_update, dash.no_update, False

        selected = _normalize_steps(run_mode, run_steps)
        if not selected:
            return dbc.Alert("No steps selected. Create run and choose steps first.", color="warning"), "", False

        sections = []
        workflow_uploaded_files = (uploaded_files or {}).copy()

        if "qc" in selected:
            sections.append(html.H5("QC"))
            sections.append(qc_callback(1, uploaded_files, project_state, samples_store, workflow_scope))
            sections.append(html.Hr())

        if "preprocess" in selected:
            # Use single explicit preprocess sample if provided
            explicit_sample_id = None
            if (uploaded_files or {}).get("preprocess_mode") == "pe":
                explicit_sample_id = preprocess_sample
            sections.extend([html.H5("Preprocess"), _run_preprocess_auto(uploaded_files, project_state, samples_store, workflow_scope, explicit_pe_sample=explicit_sample_id), html.Hr()])

        if "analysis" in selected:
            # use values from UI state (may be None)
            analysis_method = analysis_method or "vsearch_only"
            vsearch_identity = vsearch_identity or 0.97
            raxml_bootstrap = raxml_bootstrap or 100
            phylo_input_choice = phylo_input_choice or "auto"
            sections.append(html.H5("Analysis"))
            sections.append(run_analysis(1, analysis_method, vsearch_identity, uploaded_files, project_state, samples_store, workflow_scope))
            sections.append(html.Hr())

        if "analysis" in selected and "db" in selected:
            analysis_base = _step_base(project_state, "analysis")
            dada2_asv = os.path.join(analysis_base, "dada2", "repseqs.fasta")
            if not (os.path.exists(dada2_asv) and os.path.getsize(dada2_asv) > 0):
                dada2_asv = _resolve_output_path(project_state, "analysis", "final", "asv.fasta")
            if os.path.exists(dada2_asv) and os.path.getsize(dada2_asv) > 0:
                current_db_inputs = workflow_uploaded_files.get("databases", [])
                if not current_db_inputs:
                    workflow_uploaded_files["databases"] = [dada2_asv]

        if "db" in selected:
            db_outputs = [html.H5("Databases")]
            selected_dbs = db_selection or []
            blast_max_hits = blast_max_hits or 5
            blast_evalue = blast_evalue or "1e-10"
            bold_db = bold_db or 1
            bold_mode = bold_mode or 1
            if "bold" in selected_dbs:
                db_outputs.append(run_boldigger_callback(1, bold_db, bold_mode, bold_thresholds, workflow_uploaded_files, project_state))
            non_bold = [db for db in selected_dbs if db != "bold"]
            if non_bold:
                db_outputs.append(
                    run_db_classification(
                        1,
                        non_bold,
                        custom_db,
                        blast_max_hits,
                        blast_evalue,
                        workflow_uploaded_files,
                        project_state,
                        samples_store,
                        workflow_scope,
                    )
                )
            if not selected_dbs:
                db_outputs.append(html.Div("No databases selected for DB step.", style={"color": "orange"}))
            db_outputs.append(html.Hr())
            sections.extend(db_outputs)

        if "phylo" in selected:
            raxml_bootstrap = 100
            phylo_input_choice = 'auto'  # Use automatic selection in workflow mode
            sections.extend([html.H5("Phylogenetics"), run_phylo(1, raxml_bootstrap, phylo_input_choice, uploaded_files, project_state, workflow_scope), html.Hr()])

        stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        status = html.Span(f"Workflow executed at {stamp}.", style={"color": "green"})

        # PHASE 2: Apply cleanup and generate export after workflow completes
        project_dir = project_state.get("project_dir") if project_state else None
        run_dir = project_state.get("run_dir") if project_state else None
        if project_dir and run_dir:
            # Apply cleanup (old runs + temp files)
            _apply_cleanup_post_run(project_dir, run_dir)

            # Generate export ZIP
            export_zip = _generate_run_export_zip(project_state, run_config, samples_store, workflow_scope)
            if export_zip and os.path.exists(export_zip):
                zip_name = os.path.basename(export_zip)
                zip_size = os.path.getsize(export_zip) / (1024 * 1024)
                status = html.Div([
                    html.Span(f"Workflow executed at {stamp}.", style={"color": "green"}),
                    html.Br(),
                    html.Span([
                        html.A(
                            f"Export created: {zip_name}",
                            href=_download_href(export_zip),
                            download=zip_name,
                            target="_blank",
                            style={"color": "blue", "fontSize": "0.9em"},
                        ),
                        html.Span(f" ({zip_size:.1f} MB)", style={"color": "blue", "fontSize": "0.9em"})
                    ])
                ])
            elif _LAST_EXPORT_ERROR:
                status = html.Div([
                    html.Span(f"Workflow executed at {stamp}.", style={"color": "green"}),
                    html.Br(),
                    html.Span(_LAST_EXPORT_ERROR, style={"color": "orange", "fontSize": "0.9em"})
                ])
        
        return html.Div(sections), status, True
    except Exception as e:
        details = traceback.format_exc()
        return html.Div([
            dbc.Alert(f"Workflow execution error: {e}", color="danger"),
            dbc.Accordion([
                dbc.AccordionItem(html.Pre(details), title="Technical details")
            ], start_collapsed=True, always_open=False)
        ]), html.Span("Workflow failed.", style={"color": "red"}), False

@du.callback(
    output=[
        Output('upload-status-qc', 'children'),
        Output('workflow-progress', 'data', allow_duplicate=True),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-qc'
)
def get_a_list_qc(filenames):
    try:
        if filenames is None or not isinstance(filenames, (list, str)):
            return dash.no_update, dash.no_update, dash.no_update

        normalized_files = _normalize_uploaded_paths(filenames)
        if not normalized_files:
            status_children = html.Div("No files uploaded.")
            updated_uploaded_files = _merge_uploaded_files_from_ctx({'qc': []})
            return status_children, dash.no_update, updated_uploaded_files

        status_children = html.Ul([html.Li(f"Successfully uploaded: {os.path.basename(x)}") for x in normalized_files])
        progress_update = {
            'qc_done': True,
            'preprocess_done': False,
            'analysis_done': False,
            'db_done': False,
            'phylo_done': False,
            'results_done': False
        }
        uploaded_files_update = _merge_uploaded_files_from_ctx({'qc': normalized_files})
        return status_children, progress_update, uploaded_files_update
    except Exception as e:
        return html.Div(f"QC upload processing error: {e}", style={"color": "red"}), dash.no_update, dash.no_update
    return dash.no_update, dash.no_update, dash.no_update

@du.callback(
    output=[
        Output('upload-status-analysis', 'children'),
        Output('workflow-progress', 'data', allow_duplicate=True),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-analysis'
)
def get_a_list_analysis(filenames):
    try:
        if filenames is None or not isinstance(filenames, (list, str)):
            return dash.no_update, dash.no_update, dash.no_update

        normalized_files = _normalize_uploaded_paths(filenames)
        uploaded_files_update = _merge_uploaded_files_from_ctx({})
        status_message = html.Div("No files uploaded for analysis.")

        progress_update = {
            'qc_done': True,
            'preprocess_done': True,
            'analysis_done': False,
            'db_done': False,
            'phylo_done': False,
            'results_done': False
        }

        if normalized_files:
            actual_filename = normalized_files[0]
            uploaded_files_update['analysis'] = [actual_filename]
            status_message = html.Ul([html.Li(f"File for Analysis: {os.path.basename(actual_filename)}")])
            progress_update['analysis_done'] = True
            return status_message, progress_update, uploaded_files_update
        uploaded_files_update['analysis'] = []
        return status_message, dash.no_update, uploaded_files_update
    except Exception as e:
        return html.Div(f"Analysis upload processing error: {e}", style={"color": "red"}), dash.no_update, dash.no_update
    return dash.no_update, dash.no_update, dash.no_update

@du.callback(
    output=[
        Output('upload-status-databases', 'children'),
        Output('workflow-progress', 'data', allow_duplicate=True),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-databases'
)
def get_a_list_databases(filenames):
    try:
        if filenames is None or not isinstance(filenames, (list, str)):
            return dash.no_update, dash.no_update, dash.no_update

        normalized_files = _normalize_uploaded_paths(filenames)
        uploaded_files_update = _merge_uploaded_files_from_ctx({})
        status_message = html.Div("No FASTA files uploaded to BOLDigger.")

        progress_update = {
            'qc_done': True,
            'preprocess_done': True,
            'analysis_done': True,
            'db_done': False,
            'phylo_done': False,
            'results_done': False
        }

        if normalized_files:
            actual_filename = normalized_files[0]
            uploaded_files_update['databases'] = [actual_filename]
            status_message = html.Ul([html.Li(f"FASTA file for BOLDigger: {os.path.basename(actual_filename)}")])
            progress_update['db_done'] = True
            return status_message, progress_update, uploaded_files_update
        uploaded_files_update['databases'] = []
        return status_message, dash.no_update, uploaded_files_update
    except Exception as e:
        return html.Div(f"Databases upload processing error: {e}", style={"color": "red"}), dash.no_update, dash.no_update
    return dash.no_update, dash.no_update, dash.no_update

@du.callback(
    output=[
        Output('upload-status-phylo', 'children'),
        Output('workflow-progress', 'data', allow_duplicate=True),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-phylo'
)
def get_a_list_phylo(filenames):
    try:
        if filenames is None or not isinstance(filenames, (list, str)):
            return dash.no_update, dash.no_update, dash.no_update

        normalized_files = _normalize_uploaded_paths(filenames)
        uploaded_files_update = _merge_uploaded_files_from_ctx({})
        status_message = html.Div("No FASTA files uploaded for phylogenetic analysis.")

        progress_update = {
            'qc_done': True,
            'preprocess_done': True,
            'analysis_done': True,
            'db_done': True,
            'phylo_done': False,
            'results_done': False
        }
        if normalized_files:
            actual_filename = normalized_files[0]
            uploaded_files_update['phylo'] = [actual_filename]
            status_message = html.Ul([html.Li(f"FASTA file for phylogeny: {os.path.basename(actual_filename)}")])
            progress_update['phylo_done'] = True
            return status_message, progress_update, uploaded_files_update
        uploaded_files_update['phylo'] = []
        return status_message, dash.no_update, uploaded_files_update
    except Exception as e:
        return html.Div(f"Phylo upload processing error: {e}", style={"color": "red"}), dash.no_update, dash.no_update
    return dash.no_update, dash.no_update, dash.no_update

@du.callback(
    output=[
        Output('temp-se-status-store', 'data'),
        Output('uploaded-files', 'data', allow_duplicate=True)
    ],
    id='du-upload-preprocess-se' 
)
def upload_preprocess_se_to_temp_store(filenames):
    try:
        if filenames is None:
            return dash.no_update, dash.no_update

        normalized_files = _normalize_uploaded_paths(filenames)
        try:
            store_data_main = (dash.callback_context.states or {}).get('uploaded-files.data') or {}
        except Exception:
            store_data_main = {}
        if not isinstance(store_data_main, dict):
            store_data_main = {}

        updated_uploaded_files = store_data_main.copy()
        status_content_for_store = None

        
        updated_uploaded_files['preprocess_mode'] = 'se'
        
        current_filenames = normalized_files

        updated_uploaded_files['preprocess'] = current_filenames
        
        if 'preprocess_r1' in updated_uploaded_files: del updated_uploaded_files['preprocess_r1']
        if 'preprocess_r2' in updated_uploaded_files: del updated_uploaded_files['preprocess_r2']

        if current_filenames:
            status_content_for_store = [f"File Single-End: {os.path.basename(x)}" for x in current_filenames]
        else:
            status_content_for_store = ["No single-end files uploaded or uploads cancelled."]

        return status_content_for_store, updated_uploaded_files
    except Exception as e:
        error_status_content = [f"Error uploading SE: {str(e)}"]
        return error_status_content, dash.no_update

@du.callback(
    output=Output('temp-r1-path-store', 'data'), 
    id='du-upload-preprocess-r1', 
)
def upload_preprocess_r1_to_temp_store(filenames):
    try:
        normalized_files = _normalize_uploaded_paths(filenames)
        if normalized_files:
            return normalized_files[0]
        return None
    except Exception:
        return None 


@du.callback(
    output=Output('temp-r2-path-store', 'data'), 
    id='du-upload-preprocess-r2', 
)
def upload_preprocess_r2_to_temp_store(filenames):
    try:
        normalized_files = _normalize_uploaded_paths(filenames)
        if normalized_files:
            return normalized_files[0]
        return None
    except Exception:
        return None



@app.callback(
    Output('qc-results', 'children'),
    Input('run-qc-btn', 'n_clicks'),
    State('uploaded-files', 'data'),
    State('project-state', 'data'),
    State('samples-store', 'data'),
    State('workflow-scope-samples', 'data'),
    prevent_initial_call=True
)
def qc_callback(n_clicks, uploaded_files, project_state, samples_store, workflow_scope):
    if not n_clicks:
        return ""
    fastqs = uploaded_files.get('qc', []) if uploaded_files else []
    if not fastqs:
        # Use workflow-scope-samples to iterate
        collected = []
        for sample_id, sample_info in _iter_workflow_samples(samples_store, workflow_scope):
            paths = (sample_info or {}).get("paths", {})
            if paths.get("r1"):
                collected.append(paths.get("r1"))
            elif paths.get("se"):
                collected.append(paths.get("se"))
        fastqs = collected
    qc_inputs = {
        "fastqs": list(fastqs or []),
        "workflow_scope": workflow_scope or {},
    }

    def _write_qc_manifest(status, outputs=None, error="", notes=""):
        _write_step_manifest(
            project_state,
            "qc",
            _build_step_manifest("qc", status, outputs=outputs or {}, inputs=qc_inputs, error=error, notes=notes),
        )

    if not fastqs:
        _write_qc_manifest(STEP_STATUS_SKIPPED, notes="No valid FASTQ files were available for QC.")
        return html.Div("No valid FASTQ files found for QC. Please upload a non-empty FASTQ file.")

    _write_qc_manifest(STEP_STATUS_RUNNING)
    results = []
    report_outputs = {}
    failures = []
    for input_file in fastqs:
        if not _path_is_valid(input_file):
            failures.append(f"Missing or empty input FASTQ: {input_file}")
            results.append(html.Div(f"Skipping invalid FASTQ input: {os.path.basename(str(input_file))}", style={"color": "red"}))
            continue

        base_for_subdir = os.path.splitext(os.path.basename(input_file))[0] 
        output_dir_fastqc = _resolve_output_path(project_state, "qc", "intermediate", f"fastqc_results/{base_for_subdir}")
        os.makedirs(output_dir_fastqc, exist_ok=True)
        
        stdout, stderr = run_fastqc(input_file, output_dir_fastqc)

        report_paths = _fastqc_report_paths(input_file, output_dir_fastqc)
        html_valid = _path_is_valid(report_paths["html"])
        zip_valid = _path_is_valid(report_paths["zip"])
        report_outputs[os.path.basename(input_file)] = {
            "input": input_file,
            "html": report_paths["html"] if html_valid else "",
            "zip": report_paths["zip"] if zip_valid else "",
        }

        if not html_valid and not zip_valid:
            failures.append(f"FastQC did not produce a valid report for {input_file}")
            results.append(html.Div([
                html.H5(f"FastQC failed for {os.path.basename(input_file)}"),
                html.Pre(stdout),
                html.Pre(stderr),
            ], style={"color": "red"}))
            continue

        download_target = report_paths["html"] if html_valid else report_paths["zip"]
        download_label = "Download FastQC HTML Report" if html_valid else "Download FastQC archive"
        report_url_path = _download_href(download_target)

        download_link = html.A(
            f"{download_label} for {os.path.basename(input_file)}",
            href=report_url_path,
            target="_blank",
            download=os.path.basename(download_target)
        )
        report_children = [
            html.H5("Downloads"),
            download_link,
        ]
        if html_valid:
            report_children.append(dbc.Accordion([
                dbc.AccordionItem([
                    html.H6("View FastQC HTML Report"),
                    html.Iframe(
                        src=report_url_path,
                        style={"width": "100%", "height": "600px", "border": "1px solid #ccc"}
                    )
                ], title="View FastQC Report"),
                dbc.AccordionItem([
                    html.H6("FastQC stdout"),
                    html.Pre(stdout),
                    html.H6("FastQC stderr"),
                    html.Pre(stderr),
                ], title="Log FastQC"),
            ], start_collapsed=True, always_open=False, style={"marginTop": "1em"}))
        else:
            report_children.append(dbc.Accordion([
                dbc.AccordionItem([
                    html.H6("FastQC stdout"),
                    html.Pre(stdout),
                    html.H6("FastQC stderr"),
                    html.Pre(stderr),
                ], title="Log FastQC"),
            ], start_collapsed=True, always_open=False, style={"marginTop": "1em"}))
        results.append(html.Div(report_children))

    qc_status = STEP_STATUS_SUCCESS if results and not failures else (STEP_STATUS_FAILED if failures else STEP_STATUS_SKIPPED)
    _write_qc_manifest(qc_status, outputs={"reports": report_outputs}, error="; ".join(failures))
    return html.Div(results) if results else html.Div("No valid FastQC outputs were generated.", style={"color": "orange"})

@app.callback(
    Output('qc-summary-results', 'children'),
    Input('build-qc-summary-btn', 'n_clicks'),
    State('project-state', 'data'),
    prevent_initial_call=True
)
def build_qc_summary(n_clicks, project_state):
    if not n_clicks:
        return dash.no_update

    qc_manifest = _read_step_manifest(project_state, "qc", default={}) or {}
    report_entries = ((qc_manifest.get("outputs") or {}).get("reports") or {}) if isinstance(qc_manifest, dict) else {}
    entries = []

    if isinstance(report_entries, dict):
        iterable_reports = report_entries.items()
    elif isinstance(report_entries, list):
        iterable_reports = enumerate(report_entries)
    else:
        iterable_reports = []

    for sample_name, report_info in iterable_reports:
        report_info = report_info or {}
        zip_path = report_info.get("zip") if isinstance(report_info, dict) else None
        html_path = report_info.get("html") if isinstance(report_info, dict) else None
        stats = _parse_fastqc_basic_stats_from_zip(zip_path)
        if not stats and _path_is_valid(html_path):
            report_zip_guess = os.path.splitext(html_path)[0] + ".zip"
            stats = _parse_fastqc_basic_stats_from_zip(report_zip_guess)
        if not stats:
            continue
        entries.append({
            "sample": sample_name,
            "total_sequences": stats.get("Total Sequences", ""),
            "percent_gc": stats.get("%GC", ""),
            "sequence_length": stats.get("Sequence length", "")
        })

    if not entries:
        _write_step_manifest(project_state, "qc", _build_step_manifest("qc", STEP_STATUS_SKIPPED, outputs={"summary": ""}, inputs={"summary_requested": True}, notes="No validated FastQC outputs were available for summary generation."))
        return html.Div("No FastQC data files found.", style={"color": "orange"})

    df = pd.DataFrame(entries)
    results_dir = _resolve_output_path(project_state, "analysis", "final")
    out_path = os.path.join(results_dir, "qc_summary.tsv")
    df.to_csv(out_path, sep="\t", index=False)

    _write_step_manifest(
        project_state,
        "qc",
        _build_step_manifest(
            "qc",
            STEP_STATUS_SUCCESS,
            outputs={"reports": report_entries, "summary": out_path},
            inputs={"summary_requested": True},
            notes="QC summary generated from validated FastQC outputs.",
        ),
    )

    table_preview = dash_table.DataTable(
        data=df.to_dict("records"),
        columns=[{"name": i, "id": i} for i in df.columns],
        page_size=10,
        style_table={"overflowX": "auto", "marginTop": "1em"},
        style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
        style_cell={"textAlign": "left", "padding": "5px"}
    )
    link = html.A("Download QC summary", href=_download_href(out_path), target="_blank", download=os.path.basename(out_path))
    return html.Div([link, html.Hr(), table_preview])

@app.callback(
    [Output('preprocess-uploader-container', 'children'),
     Output('uploaded-files', 'data', allow_duplicate=True),
     Output('upload-status-preprocess', 'children')],
    Input('preprocess-confirm-btn', 'n_clicks'),
    [State('preprocess-mode', 'value'),
     State('preprocess-sample-selector', 'value'),
     State('samples-store', 'data'),
     State('workflow-scope-samples', 'data'),
     State('uploaded-files', 'data'),
     State('project-state', 'data')],
    prevent_initial_call=True
)
def show_preprocess_uploaders_and_reset_state(n_clicks, new_mode, preprocess_sample, samples_store, workflow_scope, current_uploaded_files_data, project_state):
    if not n_clicks: 
        return dash.no_update, dash.no_update, dash.no_update

    # Fallback if preprocess-mode component doesn't exist in DOM
    if new_mode is None:
        new_mode = (current_uploaded_files_data or {}).get('preprocess_mode', 'se')

    # Phase 1: Validate PE samples are in workflow-scope if workflow scope is set
    workflow_samples = _get_workflow_scope_samples(workflow_scope)
    if workflow_samples:  # If workflow-scope is populated
        if new_mode == "pe":
            if preprocess_sample and preprocess_sample not in workflow_samples:
                _write_preprocess_manifest(
                    STEP_STATUS_FAILED,
                    error=f"Selected preprocess sample '{preprocess_sample}' is not in workflow scope.",
                )
                return (html.Div(f"ERROR: Selected preprocess sample '{preprocess_sample}' not in workflow scope", style={"color": "red"}),
                        dash.no_update, dash.no_update)
    
    updated_uploaded_files = (current_uploaded_files_data or {}).copy()
    updated_uploaded_files['preprocess_mode'] = new_mode

    preprocess_inputs = {
        "mode": new_mode,
        "selected_sample": preprocess_sample or "",
        "workflow_scope": workflow_scope or {},
        "samples_selected": _get_workflow_scope_samples(workflow_scope),
    }

    def _write_preprocess_manifest(status, outputs=None, error="", notes=""):
        _write_step_manifest(
            project_state=project_state,
            step_name="preprocess",
            manifest=_build_step_manifest("preprocess", status, outputs=outputs or {}, inputs=preprocess_inputs, error=error, notes=notes),
        )

    # Clean up old keys
    for key in ['preprocess', 'preprocess_r1', 'preprocess_r2']:
        if key in updated_uploaded_files:
            del updated_uploaded_files[key]

    status_items = []
    target_samples = _iter_workflow_samples(samples_store, workflow_scope)
    
    if not target_samples:
        msg = "No samples selected in workflow scope. Please select samples in the Workflow Sample Scope card."
        _write_preprocess_manifest(STEP_STATUS_SKIPPED, notes=msg)
        return html.Div(msg, style={"color": "orange"}), updated_uploaded_files, html.Div(msg, style={"color": "orange"})

    if new_mode == "pe":
        sample_map = samples_store or {}

        def _resolve_pe_path_from_sample(sample_id):
            if not sample_id:
                return None, None
            sample_info = sample_map.get(sample_id, {})
            paths = sample_info.get("paths", {})
            return paths.get("r1"), paths.get("r2")

        r1_path, r2_path = None, None
        if preprocess_sample:
            r1_path, r2_path = _resolve_pe_path_from_sample(preprocess_sample)

        updated_uploaded_files['preprocess_r1'] = r1_path
        updated_uploaded_files['preprocess_r2'] = r2_path

        status_items.append(html.Li(f"File R1: {os.path.basename(r1_path) if r1_path else 'Waiting...'}"))
        status_items.append(html.Li(f"File R2: {os.path.basename(r2_path) if r2_path else 'Waiting...'}"))

        uploader_ui = html.Div()
    else: # SE mode
        se_paths = []
        for sid, sample in target_samples:
            se_path = (sample or {}).get("paths", {}).get("se")
            if se_path:
                se_paths.append(se_path)
                status_items.append(html.Li(f"{sid}: {os.path.basename(se_path)}"))
            else:
                status_items.append(html.Li(f"{sid}: missing single-end FASTQ"))
        updated_uploaded_files['preprocess'] = se_paths
        
        if len(target_samples) == len(samples_store):
            status_items.insert(0, html.Li(f"Using all samples ({len(target_samples)})."))
        else:
            status_items.insert(0, html.Li(f"Using workflow-selected samples ({len(target_samples)})."))
            
        uploader_ui = html.Div()

    helper_msg = html.Div("Using FASTQ files from project dataset (no re-upload needed).", style={"color": "#6c757d"})
    _write_preprocess_manifest(
        STEP_STATUS_PENDING,
        outputs={
            "mode": new_mode,
            "preprocess_r1": updated_uploaded_files.get("preprocess_r1", ""),
            "preprocess_r2": updated_uploaded_files.get("preprocess_r2", ""),
            "preprocess": updated_uploaded_files.get("preprocess", []),
        },
        notes="Preprocess configuration confirmed; awaiting execution.",
    )
    return uploader_ui, updated_uploaded_files, html.Ul([helper_msg] + status_items)


@app.callback(
    Output("preprocess-pe-overrides-container", "style"),
    Input("uploaded-files", "data"),
    prevent_initial_call=False
)
def toggle_preprocess_pe_overrides(uploaded_files):
    preprocess_mode = (uploaded_files or {}).get('preprocess_mode', 'se')
    if preprocess_mode == "se":
        return {"display": "none"}
    return {}


@app.callback(
    [Output('upload-status-preprocess', 'children', allow_duplicate=True),
     Output('uploaded-files', 'data', allow_duplicate=True)],
    Input('preprocess-sample-selector', 'value'),
    [State('samples-store', 'data'),
     State('uploaded-files', 'data')],
    prevent_initial_call=True
)
def update_pe_selection_status(preprocess_sample_id, samples_store, current_uploaded_files_data):
    data = (current_uploaded_files_data or {}).copy()
    if data.get('preprocess_mode') != 'pe':
        return dash.no_update, dash.no_update

    sample_map = samples_store or {}

    r1_path, r2_path = None, None
    if preprocess_sample_id:
        paths = sample_map.get(preprocess_sample_id, {}).get('paths', {})
        r1_path = paths.get('r1') or paths.get('se')
        r2_path = paths.get('r2') or paths.get('se')

    data['preprocess_r1'] = r1_path
    data['preprocess_r2'] = r2_path

    helper_msg = html.Div("Using FASTQ files from project dataset (no re-upload needed).", style={"color": "#6c757d"})
    status_items = [
        html.Li(f"File R1: {os.path.basename(r1_path) if r1_path else 'Waiting...'}"),
        html.Li(f"File R2: {os.path.basename(r2_path) if r2_path else 'Waiting...'}")
    ]
    return html.Ul([helper_msg] + status_items), data


@app.callback(
    Output('upload-status-preprocess', 'children', allow_duplicate=True),
    Input('temp-se-status-store', 'data'),
    prevent_initial_call=True
)
def update_se_status_from_temp_store(status_data_from_temp_store):
    if not status_data_from_temp_store: 
        return dash.no_update
    if isinstance(status_data_from_temp_store, list):
        return html.Ul([html.Li(x) for x in status_data_from_temp_store])
    return status_data_from_temp_store


@app.callback(
    [Output('upload-status-preprocess', 'children', allow_duplicate=True),
     Output('uploaded-files', 'data', allow_duplicate=True)],
    [Input('temp-r1-path-store', 'data'),
     Input('temp-r2-path-store', 'data')],
    [State('uploaded-files', 'data')],
    prevent_initial_call=True 
)
def consolidate_pe_uploads_and_update_status(r1_path_from_temp, r2_path_from_temp, current_uploaded_files_data):
    ctx = dash.callback_context
    if not ctx.triggered_id:
        return dash.no_update, dash.no_update

    updated_data = (current_uploaded_files_data or {}).copy()
    
    # Only act if in PE mode
    if updated_data.get('preprocess_mode') != 'pe':
        return dash.no_update, dash.no_update

    # Prioritize newly uploaded files from the temp stores
    if 'temp-r1-path-store' in ctx.triggered_ids and r1_path_from_temp:
        updated_data['preprocess_r1'] = r1_path_from_temp
    if 'temp-r2-path-store' in ctx.triggered_ids and r2_path_from_temp:
        updated_data['preprocess_r2'] = r2_path_from_temp

    # Construct status message based on the potentially updated data
    r1_status = os.path.basename(updated_data.get('preprocess_r1')) if updated_data.get('preprocess_r1') else 'Waiting...'
    r2_status = os.path.basename(updated_data.get('preprocess_r2')) if updated_data.get('preprocess_r2') else 'Waiting...'
    
    status_children_content = html.Ul([
        html.Li(f"File R1: {r1_status}"),
        html.Li(f"File R2: {r2_status}")
    ])

    return status_children_content, updated_data

@app.callback(
    Output('preprocess-results', 'children'),
    [Input('run-trim-btn', 'n_clicks'), Input('run-pear-btn', 'n_clicks')],
    State('uploaded-files', 'data'),
    State('project-state', 'data'),
    State('samples-store', 'data'),
    State('workflow-scope-samples', 'data'),
    prevent_initial_call=True
)
def preprocess_tools_callback(n_trim, n_pear, uploaded_files, project_state, samples_store, workflow_scope):
    ctx = dash.callback_context
    if not ctx.triggered:
        return ""
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    def _write_preprocess_manifest(status, outputs=None, error="", notes=""):
        _write_step_manifest(
            project_state,
            "preprocess",
            _build_step_manifest(
                "preprocess",
                status,
                outputs=outputs or {},
                inputs={
                    "mode": (uploaded_files or {}).get("preprocess_mode", "se"),
                    "workflow_scope": workflow_scope or {},
                    "samples_selected": _get_workflow_scope_samples(workflow_scope),
                },
                error=error,
                notes=notes,
            ),
        )

    if not uploaded_files:
        _write_preprocess_manifest(STEP_STATUS_SKIPPED, notes="No preprocess inputs were present in the application state.")
        return html.Div("No files uploaded in the application status.")

    mode = uploaded_files.get('preprocess_mode', 'se')

    preprocess_sources = _resolve_preprocess_sources(uploaded_files, samples_store, workflow_scope)
    r1 = preprocess_sources.get("se") if mode == 'se' else preprocess_sources.get("r1")
    r2 = preprocess_sources.get("r2") if mode == 'pe' else None

    current_fastqs_for_this_step = []
    if mode == 'se' and r1:
        current_fastqs_for_this_step.append(r1)
    elif mode == 'pe':
        if r1:
            current_fastqs_for_this_step.append(r1)
        if r2:
            current_fastqs_for_this_step.append(r2)

    valid_fastqs = [f for f in current_fastqs_for_this_step if f and (f.endswith('.fastq') or f.endswith('.fastq.gz'))]

    def _invalid_fastq_message(message, status=STEP_STATUS_FAILED):
        _write_preprocess_manifest(status, error=message)
        return html.Div(message, style={"color": "red" if status == STEP_STATUS_FAILED else "orange"})

    if button_id == "run-trim-btn":
        if mode == 'se':
            if not r1:
                return _invalid_fastq_message("Single-End Mode: No R1 file found. Please upload the file.", status=STEP_STATUS_SKIPPED)
            if r1 not in valid_fastqs:
                return _invalid_fastq_message(f"Single-End Mode: File R1 ({os.path.basename(r1) if r1 else 'N/A'}) is not a valid FASTQ.")

            input_file = r1
            base = os.path.splitext(os.path.basename(input_file))[0]
            output_file_name = f"{base}_trimmed.fastq"
            output_file_path = _resolve_output_path(project_state, "preprocess", "intermediate", output_file_name)

            stdout_trim, stderr_trim = run_trimmomatic(input_file, output_file_path, TRIMMOMATIC_JAR, TRIMMOMATIC_ADAPTERS)
            output_valid = _path_is_valid(output_file_path)
            if not output_valid:
                _write_preprocess_manifest(STEP_STATUS_FAILED, outputs={"trimmed_fastq": ""}, error="Trimmomatic did not produce a valid trimmed FASTQ.")
                return html.Div([
                    html.H5("Downloads Trimmomatic (Single-End)"),
                    html.P("Trimmed output was not created or is empty.", style={"color": "red"}),
                    dbc.Accordion([
                        dbc.AccordionItem(html.Pre(read_file_preview(output_file_path)), title="View Trimmed FASTQ"),
                        dbc.AccordionItem([html.H6("Trimmomatic stdout"), html.Pre(stdout_trim), html.H6("Trimmomatic stderr"), html.Pre(stderr_trim)], title="Log Trimmomatic"),
                    ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
                ])

            _write_preprocess_manifest(STEP_STATUS_SUCCESS, outputs={"trimmed_fastq": output_file_path})
            return html.Div([
                html.H5("Downloads Trimmomatic (Single-End)"),
                html.A("Download trimmed (single-end)", href=_download_href(output_file_path), target="_blank", download=output_file_name), 
                dbc.Accordion([
                    dbc.AccordionItem(html.Pre(read_file_preview(output_file_path)), title="View Trimmed FASTQ"),
                    dbc.AccordionItem([html.H6("Trimmomatic stdout"), html.Pre(stdout_trim), html.H6("Trimmomatic stderr"), html.Pre(stderr_trim)], title="Log Trimmomatic"),
                ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
            ])

        if mode == 'pe':
            if not r1 or not r2:
                return _invalid_fastq_message("Paired-End Mode: R1 and R2 files are required. Please upload both.", status=STEP_STATUS_SKIPPED)
            if r1 == r2:
                return _invalid_fastq_message("Paired-End Mode: R1 and R2 files cannot be the same.")
            if not (r1 in valid_fastqs and r2 in valid_fastqs):
                return _invalid_fastq_message("Paired-End Mode: One or both files are not valid FASTQ files.")

            base1 = os.path.splitext(os.path.basename(r1))[0]
            base2 = os.path.splitext(os.path.basename(r2))[0]
            out_r1_paired = _resolve_output_path(project_state, "preprocess", "intermediate", f"{base1}_trimmed_R1_paired.fastq")
            out_r1_unpaired = _resolve_output_path(project_state, "preprocess", "intermediate", f"{base1}_trimmed_R1_unpaired.fastq")
            out_r2_paired = _resolve_output_path(project_state, "preprocess", "intermediate", f"{base2}_trimmed_R2_paired.fastq")
            out_r2_unpaired = _resolve_output_path(project_state, "preprocess", "intermediate", f"{base2}_trimmed_R2_unpaired.fastq")

            stdout_trim, stderr_trim = run_trimmomatic_pe(
                r1, r2, out_r1_paired, out_r1_unpaired, out_r2_paired, out_r2_unpaired,
                TRIMMOMATIC_JAR, TRIMMOMATIC_ADAPTERS_PE
            )
            output_map = {
                "trimmed_r1_paired": out_r1_paired if _path_is_valid(out_r1_paired) else "",
                "trimmed_r1_unpaired": out_r1_unpaired if _path_is_valid(out_r1_unpaired) else "",
                "trimmed_r2_paired": out_r2_paired if _path_is_valid(out_r2_paired) else "",
                "trimmed_r2_unpaired": out_r2_unpaired if _path_is_valid(out_r2_unpaired) else "",
            }
            if not any(output_map.values()):
                _write_preprocess_manifest(STEP_STATUS_FAILED, outputs=output_map, error="Trimmomatic did not produce any valid PE outputs.")
                return html.Div([
                    html.H5("Downloads Trimmomatic (Paired-End)"),
                    html.P("No trimmed PE outputs were created.", style={"color": "red"}),
                    dbc.Accordion([
                        dbc.AccordionItem(html.Pre(read_file_preview(out_r1_paired)), title="View R1 Paired Trimmed"),
                        dbc.AccordionItem(html.Pre(read_file_preview(out_r2_paired)), title="View R2 Paired Trimmed"),
                        dbc.AccordionItem([html.H6("Trimmomatic stdout"), html.Pre(stdout_trim), html.H6("Trimmomatic stderr"), html.Pre(stderr_trim)], title="Log Trimmomatic"),
                    ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
                ])

            _write_preprocess_manifest(STEP_STATUS_SUCCESS, outputs=output_map)
            return html.Div([
                html.H5("Downloads Trimmomatic (Paired-End)"),
                html.A("Download trimmed R1 (paired)", href=_download_href(out_r1_paired), target="_blank", download=os.path.basename(out_r1_paired)), 
                html.Br(),
                html.A("Download trimmed R2 (paired)", href=_download_href(out_r2_paired), target="_blank", download=os.path.basename(out_r2_paired)), 
                dbc.Accordion([
                    dbc.AccordionItem(html.Pre(read_file_preview(out_r1_paired)), title="View R1 Paired Trimmed"),
                    dbc.AccordionItem(html.Pre(read_file_preview(out_r2_paired)), title="View R2 Paired Trimmed"),
                    dbc.AccordionItem([html.H6("Trimmomatic stdout"), html.Pre(stdout_trim), html.H6("Trimmomatic stderr"), html.Pre(stderr_trim)], title="Log Trimmomatic"),
                ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
            ])

        _write_preprocess_manifest(STEP_STATUS_FAILED, error="Unexpected preprocess mode for Trimmomatic.")
        return html.Div("Unexpected error: Invalid file configuration for Trimmomatic.")

    if button_id == "run-pear-btn":
        if mode != 'pe': 
            _write_preprocess_manifest(STEP_STATUS_SKIPPED, notes="PEAR was requested while preprocess mode is not paired-end.")
            return html.Div("PEAR can only be run in Paired-End mode.")

        if not (r1 and r2): 
            _write_preprocess_manifest(STEP_STATUS_FAILED, error="PEAR input files were not found in the current state.")
            return html.Div("PEAR: Original R1 and R2 files not found in the state. Run the upload and Trimmomatic first.")

        base1_orig = os.path.splitext(os.path.basename(r1))[0]
        base2_orig = os.path.splitext(os.path.basename(r2))[0]
        trimmed_r1_for_pear = _resolve_output_path(project_state, "preprocess", "intermediate", f"{base1_orig}_trimmed_R1_paired.fastq")
        trimmed_r2_for_pear = _resolve_output_path(project_state, "preprocess", "intermediate", f"{base2_orig}_trimmed_R2_paired.fastq")

        if not (os.path.exists(trimmed_r1_for_pear) and os.path.exists(trimmed_r2_for_pear)):
            _write_preprocess_manifest(STEP_STATUS_FAILED, error="Trimmed paired FASTQ files were not found for PEAR.")
            return html.Div(f"PEAR: Trimmed files R1 ({os.path.basename(trimmed_r1_for_pear)}) and/or R2 ({os.path.basename(trimmed_r2_for_pear)}) not found. Run Trimmomatic in paired-end mode first.")

        output_dir_pear = _resolve_output_path(project_state, "preprocess", "intermediate", f"pear_results/{base1_orig}_{base2_orig}")
        os.makedirs(output_dir_pear, exist_ok=True)

        stdout_pear, stderr_pear = run_pear(trimmed_r1_for_pear, trimmed_r2_for_pear, output_dir_pear, PEAR_PATH)

        merged_file_name = "merged.assembled.fastq"
        merged_file_path = os.path.join(output_dir_pear, merged_file_name)
        merged_valid = _path_is_valid(merged_file_path)
        output_map = {
            "merged_fastq": merged_file_path if merged_valid else "",
            "trimmed_r1_paired": trimmed_r1_for_pear if _path_is_valid(trimmed_r1_for_pear) else "",
            "trimmed_r2_paired": trimmed_r2_for_pear if _path_is_valid(trimmed_r2_for_pear) else "",
        }

        if not merged_valid:
            stderr_pear += "\n[ERROR] File merged.assembled.fastq was not created by PEAR."
            _write_preprocess_manifest(STEP_STATUS_FAILED, outputs=output_map, error="PEAR did not produce a valid merged FASTQ.")
            return html.Div([
                html.H5("Downloads PEAR"),
                html.P("Merged PEAR output was not created or is empty.", style={"color": "red"}),
                dbc.Accordion([
                    dbc.AccordionItem(html.Pre(read_file_preview(merged_file_path)), title="View Merged FASTQ"),
                    dbc.AccordionItem([html.H6("PEAR stdout"), html.Pre(stdout_pear), html.H6("PEAR stderr"), html.Pre(stderr_pear)], title="Log PEAR"),
                ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
            ])

        _write_preprocess_manifest(STEP_STATUS_SUCCESS, outputs=output_map)
        return html.Div([
            html.H5("Downloads PEAR"),
            html.A("Download merged PEAR output", href=_download_href(merged_file_path), target="_blank", download=merged_file_name), 
            dbc.Accordion([
                dbc.AccordionItem(html.Pre(read_file_preview(merged_file_path)), title="View Merged FASTQ"),
                dbc.AccordionItem([html.H6("PEAR stdout"), html.Pre(stdout_pear), html.H6("PEAR stderr"), html.Pre(stderr_pear)], title="Log PEAR"),
            ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
        ])

    _write_preprocess_manifest(STEP_STATUS_SKIPPED, notes="No preprocess action matched the triggered button.")
    return ""

@app.callback(
    Output('analysis-results', 'children'),
    Input('run-analysis-btn', 'n_clicks'),
    State('analysis-method', 'value'),
    State('vsearch-identity', 'value'),
    State('uploaded-files', 'data'),
    State('project-state', 'data'),
    State('samples-store', 'data'),
    State('workflow-scope-samples', 'data'),
    prevent_initial_call=True
)
def run_analysis(n_clicks, method, identity, uploaded_files, project_state, samples_store, workflow_scope):
    if not n_clicks or not method:
        return ""
    analysis_mode = str(method or "vsearch_only")
    methods = _analysis_mode_to_methods(analysis_mode)

    selected_ids = _get_workflow_scope_samples(workflow_scope)
    if not selected_ids:
        selected_ids = list((samples_store or {}).keys())

    manifest_rows = []
    for sid in selected_ids:
        sample = (samples_store or {}).get(sid, {})
        paths = sample.get("paths", {}) if isinstance(sample, dict) else {}
        manifest_rows.append(
            {
                "internal_sample_id": sample.get("internal_sample_id") or sample.get("pair_key") or sid,
                "display_label": sample.get("display_label") or sample.get("display_name") or sid,
                "marker": sample.get("marker"),
                "mode": sample.get("mode") or sample.get("kind") or "SE",
                "forward": paths.get("r1"),
                "reverse": paths.get("r2"),
                "single": paths.get("se"),
            }
        )

    valid_manifest, manifest_error = validate_manifest_rows(manifest_rows, analysis_mode)
    if not valid_manifest:
        _write_step_manifest(
            project_state,
            "analysis",
            _build_step_manifest(
                "analysis",
                STEP_STATUS_FAILED,
                outputs={},
                inputs={"mode": analysis_mode, "methods": methods},
                error=manifest_error,
            ),
        )
        return html.Div(manifest_error, style={"color": "red"})

    if "vsearch" in methods and identity is None:
        msg = "VSEARCH clustering cutoff (perc_identity) is required."
        _write_step_manifest(
            project_state,
            "analysis",
            _build_step_manifest("analysis", STEP_STATUS_FAILED, outputs={}, inputs={"mode": analysis_mode, "methods": methods}, error=msg),
        )
        return html.Div(msg, style={"color": "red"})
    # Create a run directory and submit background job
    project_dir = project_state.get('project_dir') if project_state else None
    run_dir = _create_run(project_dir or _get_active_root(project_state))
    description = f"analysis:{analysis_mode}:{datetime.now().isoformat()}"

    # Submit the analysis worker as a background job
    run_id = jobs.submit_job(
        target=_analysis_worker,
        args=(analysis_mode, identity, uploaded_files, project_state, samples_store, workflow_scope, run_dir),
        description=description,
    )

    # Write initial manifest file pointing to run id
    _write_step_manifest(project_state, "analysis", _build_step_manifest("analysis", STEP_STATUS_RUNNING, inputs={"mode": analysis_mode, "methods": methods, "run_id": run_id}, outputs={}))

    # Immediate UI response with run_id and basic links to job log/status
    job_status_link = f"/jobs/{run_id}/status"
    job_log_link = f"/jobs/{run_id}/log"
    return html.Div([
        html.H5("Analysis job submitted"),
        html.P(f"Run ID: {run_id}"),
        html.A("View job status", href=job_status_link, target="_blank"), html.Br(),
        html.A("Download job log", href=job_log_link, target="_blank"),
        html.P("The job is running in background. Refresh this panel to see final results."),
    ])

@app.callback(
    Output('boldigger-results', 'children'),
    Input('run-boldigger-btn', 'n_clicks'),
    State('bold-db-dropdown', 'value'),
    State('bold-mode-dropdown', 'value'),
    State('bold-thresholds', 'value'),
    State('uploaded-files', 'data'),
    State('project-state', 'data'),
    prevent_initial_call=True
)
def run_boldigger_callback(n_clicks, db_nr, mode, thresholds, uploaded_files, project_state):
    if not n_clicks:
        return ""
    
    fasta_path_input_step, _, _ = _resolve_step_input([
        ("db_upload", next((fp for fp in _normalize_uploaded_paths((uploaded_files or {}).get('databases', [])) if fp.lower().endswith(('.fasta', '.fa'))), None)),
        ("analysis_vsearch", os.path.join(_step_outputs_base(project_state, "analysis"), "vsearch", "otus.fasta")),
        ("analysis_dada2", _resolve_output_path(project_state, "analysis", "intermediate", "dada2/asv.fasta")),
    ])

    if not fasta_path_input_step or not os.path.exists(fasta_path_input_step) or os.path.getsize(fasta_path_input_step) == 0:
        return html.Div("No FASTA file found for BOLDigger. Run VSEARCH first or upload a FASTA file in the 'Reference Databases' tab.")


    db_base = _step_base(project_state, "db")
    boldigger_input_fasta_dir = os.path.join(db_base, "boldigger_input")
    os.makedirs(boldigger_input_fasta_dir, exist_ok=True)
    
    sanitized_input_basename = "boldigger_input_sanitized" 
    fasta_path_unique_sanitized = os.path.join(boldigger_input_fasta_dir, f"{sanitized_input_basename}.fasta")
    temp_unique_path = os.path.join(boldigger_input_fasta_dir, "temp_unique_for_boldigger.fasta")

    try:
        uniquify_fasta(fasta_path_input_step, temp_unique_path)
        sanitize_ids(temp_unique_path, fasta_path_unique_sanitized)
        if os.path.exists(temp_unique_path):
            os.remove(temp_unique_path)
    except Exception as e:
        return html.Div(f"Error preparing FASTA file for BOLDigger (uniquify/sanitize): {e}")

    if not os.path.exists(fasta_path_unique_sanitized) or os.path.getsize(fasta_path_unique_sanitized) == 0:
        return html.Div("Error creating sanitized FASTA file for BOLDigger. The file is empty or was not created.")

    # Cache key is based on sanitized query content + BOLD options.
    cache_key = None
    try:
        cache_key = _bold_cache_key(fasta_path_unique_sanitized, db_nr, mode, thresholds)
        cache_hit = _resolve_bold_cache_hit(cache_key)
    except Exception:
        cache_hit = None

    output_file_basename = f"{sanitized_input_basename}_identification_result"
    output_xlsx_name = f"{output_file_basename}.xlsx"
    output_parquet_name = f"{output_file_basename}.parquet"
    boldigger_outputs_dir = _step_output_path(project_state, "db", "boldigger3")
    os.makedirs(boldigger_outputs_dir, exist_ok=True)
    boldigger_logs_dir = os.path.join(boldigger_outputs_dir, "logs")
    os.makedirs(boldigger_logs_dir, exist_ok=True)
    boldigger_log_path = os.path.join(boldigger_logs_dir, "boldigger3.log")
    boldigger_manifest_path = os.path.join(boldigger_outputs_dir, "boldigger3_manifest.json")
    output_xlsx_path = os.path.join(boldigger_outputs_dir, output_xlsx_name)
    output_parquet_path = os.path.join(boldigger_outputs_dir, output_parquet_name)

    if cache_hit:
        if cache_hit.get("xlsx"):
            shutil.copy2(cache_hit["xlsx"], output_xlsx_path)
        if cache_hit.get("parquet"):
            shutil.copy2(cache_hit["parquet"], output_parquet_path)

    python_bin = _resolve_boldigger_python()
    module_ok, module_error = _check_python_module(python_bin, "boldigger3")
    if not module_ok and not cache_hit:
        manifest_payload = {
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "status": STEP_STATUS_FAILED,
            "error": f"BOLDigger3 is not available in the active Python environment ({python_bin}). {module_error}",
            "python": python_bin,
            "db": db_nr,
            "mode": mode,
            "thresholds": str(thresholds or "").strip(),
        }
        _write_json(boldigger_manifest_path, manifest_payload)
        _write_text_file(boldigger_log_path, manifest_payload["error"])
        return html.Div([
            html.H5("BOLDigger3 Output"),
            html.P("BOLDigger3 is not installed in the active Python environment.", style={"color": "red"}),
            html.P(f"Use BOLDIGGER_PYTHON or BOLDIGGER_VENV_PYTHON to point to another interpreter. Current: {python_bin}", style={"color": "#6c757d"}),
            html.Pre(module_error or "boldigger3 import check failed.")
        ])

    cmd = [
        python_bin, "-m", "boldigger3", "identify",
        "--db", str(db_nr),
        "--mode", str(mode),
        fasta_path_unique_sanitized 
    ]

    if thresholds and str(thresholds).strip():
        try:
            threshold_values = [str(int(val)) for val in str(thresholds).strip().split()]
            if threshold_values:
                cmd += ["--thresholds"] + threshold_values
        except ValueError:
            return html.Div([
                html.P("Invalid thresholds. Use numbers separated by spaces, e.g. 99 97 90 85 80", style={"color": "red"}),
                dbc.Accordion([
                    dbc.AccordionItem([
                        html.H6("BOLDigger3 stdout"),
                        html.Pre("N/A - Input error"),
                        html.H6("BOLDigger3 stderr"),
                        html.Pre("User error in configuring thresholds."),
                    ], title="Log BOLDigger"),
                ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})
            ])

    cache_used = bool(cache_hit)
    result = None
    log_accordion = None

    if not cache_used:
        boldigger_results_dir = os.path.join(db_base, "boldigger3_results")
        os.makedirs(boldigger_results_dir, exist_ok=True)

        try:
            timeout_value = int(BOLDIGGER_TIMEOUT_SEC)
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
                cwd=boldigger_results_dir,
                timeout=timeout_value,
            )
        except subprocess.TimeoutExpired:
            timeout_message = f"BOLDigger timed out after {timeout_value}s."
            _write_text_file(boldigger_log_path, timeout_message)
            _write_json(boldigger_manifest_path, {
                "created_at": datetime.now().isoformat(timespec="seconds"),
                "status": STEP_STATUS_FAILED,
                "error": timeout_message,
                "python": python_bin,
                "db": db_nr,
                "mode": mode,
                "thresholds": str(thresholds or "").strip(),
            })
            return html.Div([
                html.H5("BOLDigger3 Output"),
                html.P(f"BOLDigger timed out after {BOLDIGGER_TIMEOUT_SEC}s.", style={"color": "red"}),
                html.P("Tip: reduce input size, use cache (same query/options), or increase BOLDIGGER_TIMEOUT_SEC.", style={"color": "#6c757d"})
            ])
        except Exception as e:
            _write_text_file(boldigger_log_path, str(e))
            _write_json(boldigger_manifest_path, {
                "created_at": datetime.now().isoformat(timespec="seconds"),
                "status": STEP_STATUS_FAILED,
                "error": str(e),
                "python": python_bin,
                "db": db_nr,
                "mode": mode,
                "thresholds": str(thresholds or "").strip(),
            })
            return html.Div(f"Failed to execute BOLDigger: {str(e)}")

        log_accordion = dbc.Accordion([
            dbc.AccordionItem([
                html.H6("BOLDigger3 stdout"),
                html.Pre(result.stdout if result and result.stdout else "No stdout."),
                html.H6("BOLDigger3 stderr"),
                html.Pre(result.stderr if result and result.stderr else "No stderr."),
            ], title="Log BOLDigger"),
        ], start_collapsed=True, always_open=False, style={"marginTop": "1em"})

        if result and result.returncode != 0:
            _write_text_file(boldigger_log_path, "\n".join([result.stdout or "", result.stderr or ""]))
            _write_json(boldigger_manifest_path, {
                "created_at": datetime.now().isoformat(timespec="seconds"),
                "status": STEP_STATUS_FAILED,
                "error": f"BOLDigger3 returned code {result.returncode}",
                "python": python_bin,
                "db": db_nr,
                "mode": mode,
                "thresholds": str(thresholds or "").strip(),
            })
            if "No object named additional_data" in (result.stderr or "") or \
               "No matches found for your query sequences" in (result.stderr or ""):
                return html.Div([
                    html.H5("BOLDigger3 Output"),
                    html.P("No results were found for your sequences in BOLD. "
                           "This can happen if the sequences are too divergent, short, "
                           "or do not represent typical regions of the BOLD. "
                           "Try using less restrictive thresholds, another search mode, "
                           "or check that your sequences are suitable for identification in BOLD.", style={"color": "orange"}),
                    log_accordion
                ])
            return html.Div([
                html.H5("BOLDigger3 Output"),
                html.P(f"BOLDigger3 error (return code {result.returncode}). Check the logs below.", style={"color": "red"}),
                log_accordion
            ])
        elif not result:
            return html.Div([
                html.H5("BOLDigger3 Output"),
                html.P("BOLDigger execution failed to produce a result object.", style={"color": "red"}),
                log_accordion
            ])

    xlsx_exists = os.path.exists(output_xlsx_path) and os.path.getsize(output_xlsx_path) > 0
    parquet_exists = os.path.exists(output_parquet_path) and os.path.getsize(output_parquet_path) > 0
    
    main_content = [html.H5("BOLDigger3 Output")]
    if cache_used:
        main_content.append(html.P("Loaded from local BOLD cache (no external request needed).", style={"color": "green"}))
    else:
        main_content.append(html.P("Computed via BOLDigger and stored in local cache.", style={"color": "#0b5ed7"}))

    if not xlsx_exists and not parquet_exists:
         main_content.append(html.P(f"BOLDigger3 ran, but no main result file ({output_xlsx_name} or {output_parquet_name}) was found in the expected location: {boldigger_input_fasta_dir}. Check the logs.", style={"color": "orange"}))
    else:
        if xlsx_exists:
            main_content.append(html.A(f"Download BOLDigger3 Results ({output_xlsx_name})", href=_download_href(output_xlsx_path), download=output_xlsx_name, target="_blank"))
            main_content.append(html.Br())
            
            try:
                df = pd.read_excel(output_xlsx_path)
                table_preview = dash_table.DataTable(
                    data=df.to_dict('records'),
                    columns=[{'name': i, 'id': i} for i in df.columns],
                    page_size=10,
                    style_table={'overflowX': 'auto', 'marginTop': '1em'},
                    style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
                    style_cell={
                        'height': 'auto',
                        'minWidth': '100px', 'width': '150px', 'maxWidth': '250px',
                        'whiteSpace': 'normal',
                        'textAlign': 'left',
                        'padding': '5px'
                    },
                    filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                )
                main_content.extend([html.Hr(), html.H6("Results Table Preview:"), table_preview])
            except Exception as e:
                main_content.append(html.P(f"Error reading or displaying Excel file: {str(e)}", style={"color": "orange"}))
        else:
            main_content.append(html.P(f"File {output_xlsx_name} not found or is empty.", style={"color": "red"}))
        
        if parquet_exists:
            main_content.append(html.A(f"Download BOLDigger3 Results ({output_parquet_name})", href=_download_href(output_parquet_path), download=output_parquet_name, target="_blank"))
            main_content.append(html.Br())
        elif not xlsx_exists: 
            main_content.append(html.P(f"File {output_parquet_name} not found or is empty.", style={"color": "red"}))

    if not cache_used and (xlsx_exists or parquet_exists) and cache_key:
        try:
            _save_bold_cache(
                cache_key,
                source_xlsx=output_xlsx_path if xlsx_exists else None,
                source_parquet=output_parquet_path if parquet_exists else None,
                metadata={
                    "created_at": datetime.now().isoformat(timespec="seconds"),
                    "db": db_nr,
                    "mode": mode,
                    "thresholds": str(thresholds or "").strip(),
                }
            )
            _write_json(boldigger_manifest_path, {
                "created_at": datetime.now().isoformat(timespec="seconds"),
                "status": STEP_STATUS_SUCCESS,
                "python": python_bin,
                "db": db_nr,
                "mode": mode,
                "thresholds": str(thresholds or "").strip(),
                "outputs": {
                    "xlsx": output_xlsx_path if xlsx_exists else "",
                    "parquet": output_parquet_path if parquet_exists else "",
                    "log": boldigger_log_path,
                },
            })
        except Exception:
            pass

    if log_accordion:
        return html.Div(main_content + [log_accordion])
    return html.Div(main_content)


@app.callback(
    Output('phylo-results', 'children'),
    Input('run-phylo-btn', 'n_clicks'),
    State('raxml-bootstrap-input', 'value'),
    State('phylo-input-choice', 'value'),
    State('uploaded-files', 'data'),
    State('project-state', 'data'),
    State('workflow-scope-samples', 'data'),
    prevent_initial_call=True
)
def run_phylo(n_clicks, n_bootstrap, phylo_input_choice, uploaded_files, project_state, workflow_scope):
    return _run_phylo_impl(n_clicks, n_bootstrap, uploaded_files, project_state, workflow_scope, phylo_input_choice)


def _resolve_analysis_output(project_state, prefer='auto'):
    manifest_path = _resolve_output_path(project_state, "analysis", "final", "analysis_manifest.json")
    def _first_valid(candidates):
        for candidate in candidates:
            if candidate and os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                return os.path.abspath(candidate)
        return None

    otu_candidates = [
        _resolve_output_path(project_state, "analysis", "final", "otus.fasta"),
        _step_output_path(project_state, "analysis", "otus.fasta"),
        os.path.join(_step_outputs_base(project_state, "analysis"), "vsearch", "otus.fasta"),
    ]
    asv_candidates = [
        _resolve_output_path(project_state, "analysis", "final", "asv.fasta"),
        _step_output_path(project_state, "analysis", "asv.fasta"),
        os.path.join(_step_outputs_base(project_state, "analysis"), "dada2", "repseqs.fasta"),
    ]

    if not os.path.exists(manifest_path):
        if prefer in ('auto', 'otu'):
            resolved = _first_valid(otu_candidates)
            if resolved:
                return resolved
        if prefer in ('auto', 'asv'):
            resolved = _first_valid(asv_candidates)
            if resolved:
                return resolved
        return None
    try:
        m = _read_json(manifest_path)
        if not m:
            return None
        methods = m.get('methods', {}) if isinstance(m, dict) else {}
        if prefer == 'otu':
            manifest_otu = _manifest_path_value(methods, 'vsearch', 'outputs', 'otus_fasta') or _manifest_path_value(methods, 'vsearch', 'otus_fasta')
            resolved = _first_valid([manifest_otu] + otu_candidates)
            if resolved:
                return resolved
        if prefer == 'asv':
            manifest_asv = _manifest_path_value(methods, 'dada2', 'outputs', 'asv_fasta') or _manifest_path_value(methods, 'dada2', 'asv_fasta')
            resolved = _first_valid([manifest_asv] + asv_candidates)
            if resolved:
                return resolved
        # auto: prefer vsearch then dada2
        manifest_otu = _manifest_path_value(methods, 'vsearch', 'outputs', 'otus_fasta') or _manifest_path_value(methods, 'vsearch', 'otus_fasta')
        resolved = _first_valid([manifest_otu] + otu_candidates)
        if resolved:
            return resolved
        manifest_asv = _manifest_path_value(methods, 'dada2', 'outputs', 'asv_fasta') or _manifest_path_value(methods, 'dada2', 'asv_fasta')
        resolved = _first_valid([manifest_asv] + asv_candidates)
        if resolved:
            return resolved
    except Exception:
        return None
    return None


def _is_valid_newick(path):
    """Returns True if file exists, is non-empty, and starts with '(' (valid Newick)."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            content = f.read(200).strip()
        return content.startswith("(")
    except Exception:
        return False


def _run_phylo_impl(n_clicks, n_bootstrap, uploaded_files, project_state, workflow_scope, input_choice):
    if not n_clicks:
        return ""
    if n_bootstrap is None or n_bootstrap < 10:
        n_bootstrap = 100 
    phylo_inputs = {
        "input_choice": input_choice or "auto",
        "n_bootstrap": n_bootstrap,
    }

    def _write_phylo_manifest(status, outputs=None, error="", notes=""):
        _write_step_manifest(
            project_state,
            "phylo",
            _build_step_manifest("phylo", status, outputs=outputs or {}, inputs=phylo_inputs, error=error, notes=notes),
        )

    _write_phylo_manifest(STEP_STATUS_RUNNING)
    # Resolve input via manifest or uploaded file
    if input_choice in ('otu', 'asv'):
        otus_fasta_input = _resolve_analysis_output(project_state, prefer=input_choice)
    else:
        otus_fasta_input = _resolve_analysis_output(project_state, prefer='auto')
        if not otus_fasta_input:
            otus_fasta_input = next((fp for fp in _normalize_uploaded_paths((uploaded_files or {}).get('phylo', [])) if fp.lower().endswith(('.fasta', '.fa'))), None)

    if not otus_fasta_input or not os.path.exists(otus_fasta_input) or os.path.getsize(otus_fasta_input) == 0:
        _write_phylo_manifest(STEP_STATUS_FAILED, error="FASTA input for phylogeny was not found or was empty.")
        return html.Div("FASTA input for phylogeny not found. Run VSEARCH or DADA2 first, or upload a FASTA file in the 'Phylogenetic Analysis' tab.")

    phylo_intermediate_dir = _resolve_output_path(project_state, "phylo", "intermediate", "phylo_intermediate")
    # ensure intermediate dir exists before writing
    try:
        os.makedirs(phylo_intermediate_dir, exist_ok=True)
    except Exception:
        pass

    otus_fasta_unique = os.path.join(phylo_intermediate_dir, "otus_unique_for_phylo.fasta")
    if not os.path.exists(otus_fasta_input) or os.path.getsize(otus_fasta_input) == 0:
        _write_phylo_manifest(STEP_STATUS_FAILED, error="Resolved phylo FASTA input became unavailable before processing.")
        return html.Div("FASTA input for phylogeny not found or empty. Run analysis to produce OTU/ASV FASTA first.")
    try:
        uniquify_fasta(otus_fasta_input, otus_fasta_unique)
    except Exception as e:
        _write_phylo_manifest(STEP_STATUS_FAILED, error=f"Error preparing FASTA for phylogeny: {e}")
        return html.Div(f"Error preparing FASTA for phylogeny: {e}", style={"color": "red"})

    aligned_fasta = os.path.join(phylo_intermediate_dir, "aligned.fasta")
    filtered_aligned_fasta = os.path.join(phylo_intermediate_dir, "aligned_filtered.fasta")
    aligned_safe_ids = os.path.join(phylo_intermediate_dir, "aligned_safe_ids.fasta")
    
    tree_file = os.path.join(phylo_intermediate_dir, "tree.nwk") 
    raxml_output_dir = _resolve_output_path(project_state, "phylo", "final", "raxmlng_results")

    mafft_ok, mafft_log = run_mafft(otus_fasta_unique, aligned_fasta)
    if not mafft_ok:
        _write_phylo_manifest(STEP_STATUS_FAILED, error=f"MAFFT alignment failed: {mafft_log}")
        return html.Div(f"MAFFT alignment failed: {mafft_log}")

    filter_empty_aligned_fasta(aligned_fasta, filtered_aligned_fasta, min_real=10)
    sanitize_ids(filtered_aligned_fasta, aligned_safe_ids)

    if not os.path.exists(aligned_safe_ids) or os.path.getsize(aligned_safe_ids) == 0:
        _write_phylo_manifest(STEP_STATUS_FAILED, error="No valid sequences remained after alignment and filtering.")
        return html.Div("No valid sequences remained after alignment and filtering for phylogenetic analysis.")

    fasttree_ok, fasttree_log = run_fasttree(aligned_safe_ids, tree_file)
    raxml_ok, raxml_log = run_raxml(aligned_safe_ids, raxml_output_dir, n_bootstrap)

    phylo_outputs_dir = _step_outputs_base(project_state, "phylo")
    canonical_aligned_safe_ids = os.path.join(phylo_outputs_dir, "aligned_safe_ids.fasta")
    canonical_fasttree_tree = os.path.join(phylo_outputs_dir, "tree.nwk")
    canonical_raxml_best_tree = os.path.join(phylo_outputs_dir, "raxml.raxml.bestTree")
    if os.path.exists(aligned_safe_ids) and os.path.getsize(aligned_safe_ids) > 0:
        shutil.copy2(aligned_safe_ids, canonical_aligned_safe_ids)
    if os.path.exists(tree_file) and os.path.getsize(tree_file) > 0:
        shutil.copy2(tree_file, canonical_fasttree_tree)
    
    raxml_best_tree_leaf_filename = "raxml.raxml.bestTree"
    raxml_best_tree_file = os.path.join(raxml_output_dir, raxml_best_tree_leaf_filename) 

    rel_raxml_best_tree_path = None
    if raxml_ok and os.path.exists(raxml_best_tree_file) and os.path.getsize(raxml_best_tree_file) > 0:
        rel_raxml_best_tree_path = raxml_best_tree_file
        shutil.copy2(raxml_best_tree_file, canonical_raxml_best_tree)

    alias_aligned_safe_ids = _resolve_output_path(project_state, "phylo", "final", "aligned_safe_ids.fasta")
    alias_fasttree_tree = _resolve_output_path(project_state, "phylo", "final", "tree.nwk")
    alias_raxml_best_tree = _resolve_output_path(project_state, "phylo", "final", "raxml.raxml.bestTree")
    _mirror_output_to_results(canonical_aligned_safe_ids, alias_aligned_safe_ids)
    _mirror_output_to_results(canonical_fasttree_tree, alias_fasttree_tree)
    _mirror_output_to_results(canonical_raxml_best_tree, alias_raxml_best_tree)

    fasttree_viz_content = html.P("FastTree failed or tree file is empty/invalid.", style={"color": "red"})
    if fasttree_ok and os.path.exists(tree_file) and os.path.getsize(tree_file) > 0:
        try:
            fasttree_obj = Phylo.read(tree_file, "newick")
            fasttree_elements = biopython_to_cytoscape_elements(fasttree_obj)
            
            root_node_id_ft = None
            if fasttree_elements and fasttree_elements[0] and 'data' in fasttree_elements[0] and 'id' in fasttree_elements[0]['data']:
                 root_node_id_ft = f"#{fasttree_elements[0]['data']['id']}"

            if fasttree_elements:
                fasttree_viz_content = cyto.Cytoscape(
                    id='cytoscape-fasttree',
                    elements=fasttree_elements,
                    layout={
                        'name': 'breadthfirst',
                        'directed': True,
                        'spacingFactor': 1.25, 
                        'grid': False,
                        'roots': root_node_id_ft
                    },
                    style={'width': '100%', 'height': '800px', 'border': '1px solid #ccc'}, 
                    stylesheet=[
                        {'selector': 'node', 'style': {'label': 'data(label)', 'font-size': '9px', 'text-valign': 'center', 'text-halign': 'center', 'background-color': '#ADD8E6', 'width': 'mapData(label.length, 0, 20, 15, 80)', 'height': 'mapData(label.length, 0, 20, 15, 80)', 'text-wrap': 'wrap', 'text-max-width': '70px'}},
                        {'selector': 'edge', 'style': {'width': 1, 'line-color': '#9dbaea', 'curve-style': 'bezier'}}, 
                        {'selector': '[label = ""]', 'style':{'width': '5px', 'height': '5px', 'background-color': 'grey', 'font-size': '7px'}}, 
                        {'selector': 'node[!label]', 'style':{'width': '5px', 'height': '5px', 'background-color': 'grey', 'font-size': '7px'}}, 
                        {'selector': 'node[label^="."]', 'style':{'font-size': '7px', 'width': '5px', 'height': '5px'}}, 
                    ]
                )
            else:
                fasttree_viz_content = html.P("Could not generate FastTree visualization elements.", style={"color": "orange"})
        except Exception as e:
            fasttree_viz_content = html.P(f"Error parsing FastTree Newick for visualization: {e}", style={"color": "red"})
    elif not fasttree_ok:
         fasttree_viz_content = html.P(f"FastTree failed: {fasttree_log}", style={"color": "red"})

    raxml_viz_content = html.P("RAxML-NG failed or tree file is empty/invalid.", style={"color": "red"})
    if raxml_ok and rel_raxml_best_tree_path and os.path.exists(raxml_best_tree_file) and os.path.getsize(raxml_best_tree_file) > 0:
        try:
            raxml_tree_obj = Phylo.read(raxml_best_tree_file, "newick")
            raxml_elements = biopython_to_cytoscape_elements(raxml_tree_obj)

            root_node_id_raxml = None
            if raxml_elements and raxml_elements[0] and 'data' in raxml_elements[0] and 'id' in raxml_elements[0]['data']:
                 root_node_id_raxml = f"#{raxml_elements[0]['data']['id']}"

            if raxml_elements:
                raxml_viz_content = cyto.Cytoscape(
                    id='cytoscape-raxml',
                    elements=raxml_elements,
                    layout={
                        'name': 'breadthfirst',
                        'directed': True,
                        'spacingFactor': 1.25, 
                        'grid': False,
                        'roots': root_node_id_raxml
                    },
                    style={'width': '100%', 'height': '800px', 'border': '1px solid #ccc'}, 
                    stylesheet=[ 
                        {'selector': 'node', 'style': {'label': 'data(label)', 'font-size': '9px', 'text-valign': 'center', 'text-halign': 'center', 'background-color': '#90EE90', 'width': 'mapData(label.length, 0, 20, 15, 80)', 'height': 'mapData(label.length, 0, 20, 15, 80)', 'text-wrap': 'wrap', 'text-max-width': '70px'}},
                        {'selector': 'edge', 'style': {'width': 1, 'line-color': '#9dbaea', 'curve-style': 'bezier'}},
                        {'selector': '[label = ""]', 'style':{'width': '5px', 'height': '5px', 'background-color': 'grey', 'font-size': '7px'}},
                        {'selector': 'node[!label]', 'style':{'width': '5px', 'height': '5px', 'background-color': 'grey', 'font-size': '7px'}},
                        {'selector': 'node[label^="."]', 'style':{'font-size': '7px', 'width': '5px', 'height': '5px'}},
                    ]
                )
            else:
                raxml_viz_content = html.P("Could not generate RAxML-NG visualization elements.", style={"color": "orange"})
        except Exception as e:
            raxml_viz_content = html.P(f"Error parsing RAxML-NG Newick for visualization: {e}", style={"color": "red"})
    elif not raxml_ok:
        raxml_viz_content = html.P(f"RAxML-NG failed. Check log.", style={"color": "red"})

    final_outputs = {
        "input_fasta": otus_fasta_input if _path_is_valid(otus_fasta_input) else "",
        "unique_fasta": otus_fasta_unique if _path_is_valid(otus_fasta_unique) else "",
        "aligned_fasta": aligned_fasta if _path_is_valid(aligned_fasta) else "",
        "filtered_aligned_fasta": filtered_aligned_fasta if _path_is_valid(filtered_aligned_fasta) else "",
        "aligned_safe_ids": canonical_aligned_safe_ids if _path_is_valid(canonical_aligned_safe_ids) else "",
        "fasttree_tree": canonical_fasttree_tree if _path_is_valid(canonical_fasttree_tree) else "",
        "raxml_best_tree": canonical_raxml_best_tree if _path_is_valid(canonical_raxml_best_tree) else "",
        "aligned_safe_ids_result": alias_aligned_safe_ids if _path_is_valid(alias_aligned_safe_ids) else "",
        "fasttree_tree_result": alias_fasttree_tree if _path_is_valid(alias_fasttree_tree) else "",
        "raxml_best_tree_result": alias_raxml_best_tree if _path_is_valid(alias_raxml_best_tree) else "",
    }
    # Stricter phylogenetic success criteria: require a valid Newick tree from FastTree or RAxML-NG
    phylo_notes = []
    fasttree_valid = _is_valid_newick(canonical_fasttree_tree)
    raxml_valid = _is_valid_newick(canonical_raxml_best_tree)
    if fasttree_valid or raxml_valid:
        phylo_status = STEP_STATUS_SUCCESS
    elif _path_is_valid(canonical_aligned_safe_ids):
        phylo_status = STEP_STATUS_FAILED
        phylo_notes.append("Alignment produced but no valid tree was generated by FastTree or RAxML-NG")
    else:
        phylo_status = STEP_STATUS_FAILED
    if not fasttree_ok:
        phylo_notes.append("FastTree failed")
    if not raxml_ok:
        phylo_notes.append("RAxML-NG failed")
    _write_phylo_manifest(phylo_status, outputs=final_outputs, notes="; ".join(phylo_notes))

    return html.Div([
        html.H5("Downloads Phylogenetics"),
        html.A("Download Aligned FASTA (Sanitized)", href=_download_href(canonical_aligned_safe_ids), download=os.path.basename(canonical_aligned_safe_ids), target="_blank"),
        html.Br(),
        html.A("Download Newick Tree (FastTree)", href=_download_href(canonical_fasttree_tree), download=os.path.basename(canonical_fasttree_tree), target="_blank") if fasttree_ok and os.path.exists(canonical_fasttree_tree) else html.P("FastTree Newick file not available.", style={"color":"orange"}),
        html.Br(),
        (html.A("Download RAxML Best Tree", href=_download_href(canonical_raxml_best_tree), download=os.path.basename(canonical_raxml_best_tree), target="_blank")
         if _path_is_valid(canonical_raxml_best_tree) else html.P("RAxML best tree not found or RAxML failed.", style={"color": "red" if not raxml_ok else "orange"})),
        html.Br(),
        html.Hr(),
        dbc.Accordion([
            dbc.AccordionItem(html.Pre(read_file_preview(canonical_aligned_safe_ids)), title="View FASTA Alignment (Sanitized)"),
            dbc.AccordionItem(fasttree_viz_content, title="View Phylogenetic Tree (FastTree)"),
            dbc.AccordionItem(raxml_viz_content, title="View Phylogenetic Tree (RAxML-NG)"),
            dbc.AccordionItem([html.H6("MAFFT Alignment Log"), html.Pre(mafft_log if mafft_log else "No MAFFT log output.")], title="Log MAFFT"),
            dbc.AccordionItem([html.H6("FastTree Log"), html.Pre(fasttree_log if fasttree_log else "No FastTree log output.")], title="Log FastTree"),
            dbc.AccordionItem([html.H6("RAxML-NG Output"), html.Pre(raxml_log if raxml_log else "No RAxML-NG log output.")], title="Log RAxML-NG"),
        ], start_collapsed=True, always_open=False, style={"marginTop": "2em"})
    ])


@app.server.route('/uploads/<path:path>')
def serve_file(path):
    return send_from_directory(UPLOAD_DIRECTORY, path)

@app.server.route('/files/<path:path>')
def serve_project_file(path):
    abs_path = os.path.abspath(os.path.join(os.getcwd(), path))
    projects_root = os.path.abspath(PROJECTS_ROOT)
    if not abs_path.startswith(projects_root + os.sep):
        return ("Not Found", 404)
    return send_from_directory(os.path.dirname(abs_path), os.path.basename(abs_path))

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('run-qc-btn', 'n_clicks'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_qc_done(n, progress):
    if n:
        progress['qc_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    [Input('run-trim-btn', 'n_clicks'), Input('run-pear-btn', 'n_clicks')],
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_preprocess_done(n_trim, n_pear, progress):
    if n_trim or n_pear:
        progress['preprocess_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('run-analysis-btn', 'n_clicks'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_analysis_done(n, progress):
    if n:
        progress['analysis_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('db-selection', 'value'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_db_done(selected, progress):
    if selected:
        progress['db_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('run-phylo-btn', 'n_clicks'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_phylo_done(n, progress):
    if n:
        progress['phylo_done'] = True
    return progress

@app.callback(
    Output('workflow-progress', 'data', allow_duplicate=True),
    Input('final-results', 'children'),
    State('workflow-progress', 'data'),
    prevent_initial_call=True
)
def mark_results_done(children, progress):
    if children:
        progress['results_done'] = True
    return progress

@app.callback(
    Output('db-info', 'children'),
    Input('db-selection', 'value'),
    State('samples-store', 'data'),
    State('workflow-scope-samples', 'data')
)
def show_db_info(selected, samples_store, workflow_scope):
    selected = selected or []
    if not selected:
        return html.Div("No databases selected.")

    marker = _infer_marker_from_samples(samples_store, workflow_scope)

    advice = []
    if marker == 'COI':
        advice.append(html.H6("Marker inferred: COI (mitochondrial COX1)"))
        advice.append(html.Ul([
            html.Li("Recommended: BOLD (primary)"),
            html.Li("Fallback: NCBI GenBank"),
            html.Li("Note: SILVA/PR2 are not appropriate for COI")
        ]))
    elif marker == '16S':
        advice.append(html.H6("Marker inferred: 16S (prokaryotic SSU)"))
        advice.append(html.Ul([
            html.Li("Recommended: SILVA (primary)"),
            html.Li("Fallback: NCBI GenBank")
        ]))
    elif marker == '18S':
        advice.append(html.H6("Marker inferred: 18S (eukaryotic SSU)"))
        advice.append(html.Ul([
            html.Li("Recommended: PR2 (primary)"),
            html.Li("Alternative: SILVA"),
            html.Li("Fallback: NCBI GenBank")
        ]))
    elif marker == 'ITS':
        advice.append(html.H6("Marker inferred: ITS"))
        advice.append(html.Ul([
            html.Li("Note: No dedicated ITS DB is bundled in this project by default."),
            html.Li("Recommended: use custom DB or NCBI as fallback")
        ]))
    else:
        advice.append(html.H6("Marker could not be confidently inferred"))
        advice.append(html.P("No clear marker detected from sample names. Showing general recommendations."))
        advice.append(html.Ul([
            html.Li("SILVA: good for rRNA/SSU"),
            html.Li("PR2: good for eukaryotic 18S"),
            html.Li("BOLD: COI barcode sequences (animals/plants)"),
            html.Li("NCBI: general fallback for many markers")
        ]))

    # Also show which of the selected DBs are present/available as quick status
    db_status_lines = []
    db_descriptions = {
        "ncbi": "NCBI GenBank (hybrid)",
        "silva": "SILVA SSU Ref NR99",
        "pr2": "PR2 5.1.1",
        "bold": "BOLD (BOLDigger3)",
        "custom": "Custom DB (project upload)",
    }
    # (db availability status can be added here if needed)

    children = [html.Div(advice)]
    # show selected DB list
    children.append(html.H6("Selected DBs"))
    children.append(html.Ul([html.Li(db_descriptions.get(db, db)) for db in selected]))
    return html.Div(children)

@app.callback(
    Output("boldigger-options-container", "style"),
    Input("db-selection", "value"),
    prevent_initial_call=False
)
def toggle_boldigger_options(selected):
    selected = selected or []
    if "bold" in selected:
        return {}
    return {"display": "none"}

@app.callback(
    Output("db-classify-results", "children"),
    Input("run-db-classify-btn", "n_clicks"),
    State("db-selection", "value"),
    State("custom-db-dropdown", "value"),
    State("blast-max-hits", "value"),
    State("blast-evalue", "value"),
    State("uploaded-files", "data"),
    State("project-state", "data"),
    State("samples-store", "data"),
    State("workflow-scope-samples", "data"),
    prevent_initial_call=True
)
def run_db_classification(n_clicks, selected_dbs, custom_db_path, max_hits, evalue, uploaded_files, project_state, samples_store, workflow_scope):
    if not n_clicks:
        return ""

    def _write_db_manifest(status, outputs=None, error="", notes=""):
        _write_step_manifest(
            project_state,
            "db",
            _build_step_manifest(
                "db",
                status,
                outputs=outputs or {},
                inputs={
                    "selected_dbs": list(selected_dbs or []),
                    "max_hits": max_hits or 5,
                    "evalue": evalue or "1e-10",
                    "custom_db_path": custom_db_path or "",
                    "workflow_scope": workflow_scope or {},
                },
                error=error,
                notes=notes,
            ),
        )

    if not selected_dbs:
        _write_db_manifest(STEP_STATUS_SKIPPED, notes="No databases were selected for classification.")
        return html.Div("No databases selected.")

    if not shutil.which("blastn"):
        _write_db_manifest(STEP_STATUS_FAILED, error="BLAST+ is not available in the current environment.")
        return html.Div("BLAST+ not found. Install ncbi-blast+ or use Docker.")

    blast_version = _blast_version_info()
    blast_logs_dir = _step_output_path(project_state, "db", "blast_logs")
    os.makedirs(blast_logs_dir, exist_ok=True)

    def _write_blast_log(db_name, stdout, stderr, attempts=None, extra_notes=None):
        log_path = os.path.join(blast_logs_dir, f"blast_{db_name}.log")
        payload = [
            f"# BLAST version: {blast_version}",
            f"# Database: {db_name}",
        ]
        if extra_notes:
            payload.append(f"# Notes: {extra_notes}")
        if attempts:
            payload.append("# Attempts:")
            payload.append(json.dumps(attempts, indent=2, ensure_ascii=False))
        payload.extend([
            "",
            "# stdout",
            stdout or "",
            "",
            "# stderr",
            stderr or "",
        ])
        _write_text_file(log_path, "\n".join(payload) + "\n")
        return log_path

    # Priority: OTUs/ASVs from analysis, then uploaded FASTA, then pooled
    query_fasta = None
    analysis_otus = _resolve_output_path(project_state, "analysis", "final", "otus.fasta")
    analysis_asv = _resolve_output_path(project_state, "analysis", "final", "asv.fasta")
    if os.path.exists(analysis_otus) and os.path.getsize(analysis_otus) > 0:
        query_fasta = analysis_otus
    elif os.path.exists(analysis_asv) and os.path.getsize(analysis_asv) > 0:
        query_fasta = analysis_asv
    else:
        # Fallback to uploaded or pooled
        upload_candidates = _normalize_uploaded_paths((uploaded_files or {}).get("databases", []))
        if upload_candidates:
            query_fasta = upload_candidates[0]
        else:
            # Build pooled from workflow scope
            pooled_dir = _resolve_output_path(project_state, "analysis", "intermediate", "pooled_analysis")
            pooled_fasta, _ = _build_pooled_fasta_for_analysis(samples_store, workflow_scope, pooled_dir)
            if pooled_fasta and os.path.exists(pooled_fasta):
                query_fasta = pooled_fasta

    # Write a small debug log for query_fasta resolution instead of printing
    _db_debug_log = os.path.join(_step_output_path(project_state, "db", "blast_logs"), "query_resolution.log")
    os.makedirs(os.path.dirname(_db_debug_log), exist_ok=True)
    _write_text_file(_db_debug_log, "\n".join([
        f"[{datetime.now().isoformat()}] query_fasta resolution",
        f"upload_candidates: {upload_candidates}",
        f"resolved query_fasta: {query_fasta}",
        f"valid: {_path_is_valid(query_fasta) if query_fasta else False}",
    ]) + "\n")

    if not query_fasta or not os.path.exists(query_fasta):
        _write_db_manifest(STEP_STATUS_SKIPPED, notes="No query FASTA was available for database classification.")
        return html.Div(
            "No final OTU/ASV FASTA found for classification. Run the analysis step first, or upload a FASTA file explicitly in the Databases step.",
            style={"color": "red"}
        )

    shared_db_root = DB_ROOT
    os.makedirs(shared_db_root, exist_ok=True)
    results_dir = os.path.join(_get_active_root(project_state), "results", "db_classification")
    os.makedirs(results_dir, exist_ok=True)
    db_outputs_dir = _step_output_path(project_state, "db", "db_classification")
    os.makedirs(db_outputs_dir, exist_ok=True)

    _update_run_session_metadata(project_state, {
        "selected_databases": list(selected_dbs),
        "db_parameters": {
            "max_hits": max_hits or 5,
            "evalue": evalue or "1e-10",
            "custom_db_path": custom_db_path or "",
        },
    })

    max_hits = max_hits or 5
    evalue = evalue or "1e-10"

    accord_items = []
    db_method_info = {}
    db_outputs = {}
    db_errors = []
    for db_name in selected_dbs:
        db_label = db_name.upper()
        output_path = os.path.join(db_outputs_dir, f"blast_{db_name}.tsv")
        attempts = []

        if db_name == "ncbi":
            if NCBI_BLAST_DB:
                stdout, stderr, code, attempts = _run_blastn_with_retry(query_fasta, NCBI_BLAST_DB, output_path, max_hits=max_hits, evalue=evalue, remote=False)
                db_method_info[db_name] = f"Source: local NCBI subset/index ({NCBI_BLAST_DB}); BLAST {blast_version}"
            elif NCBI_REMOTE_FALLBACK:
                stdout, stderr, code, attempts = _run_blastn_with_retry(query_fasta, "nt", output_path, max_hits=max_hits, evalue=evalue, remote=True)
                db_method_info[db_name] = f"Source: NCBI remote (hybrid fallback when no local subset is configured); BLAST {blast_version}"
            else:
                accord_items.append(dbc.AccordionItem(
                    html.P(
                        "NCBI selected but no local NCBI_BLAST_DB is configured and remote fallback is disabled (NCBI_REMOTE_FALLBACK=0).",
                        style={"color": "red"}
                    ),
                    title=f"{db_label} (configuration required)"
                ))
                db_errors.append("NCBI selected but no local NCBI_BLAST_DB is configured and remote fallback is disabled.")
                continue
        elif db_name == "silva":
            silva_dir = os.path.join(shared_db_root, "silva")
            fasta_path, err = ensure_fasta_db("silva", DEFAULT_SILVA_URL, silva_dir, local_path=DEFAULT_SILVA_PATH)
            if err:
                accord_items.append(dbc.AccordionItem(html.P(err, style={"color": "red"}), title=f"{db_label} (error)"))
                continue
            db_prefix = os.path.join(silva_dir, "silva")
            ok, log = ensure_blast_db(fasta_path, db_prefix)
            if not ok:
                accord_items.append(dbc.AccordionItem(html.Pre(log), title=f"{db_label} (makeblastdb failed)"))
                _write_blast_log(db_name, "", log, extra_notes=f"makeblastdb failed for {db_label}")
                db_errors.append(f"{db_label}: makeblastdb failed")
                continue
            stdout, stderr, code, attempts = _run_blastn_with_retry(query_fasta, db_prefix, output_path, max_hits=max_hits, evalue=evalue, remote=False)
        elif db_name == "pr2":
            pr2_dir = os.path.join(shared_db_root, "pr2")
            fasta_path, taxmap_path, err = prepare_pr2_mothur_db(pr2_dir)
            if (not fasta_path or not os.path.exists(fasta_path)) and not err:
                fasta_path, err = ensure_fasta_db("pr2", DEFAULT_PR2_URL, pr2_dir, local_path=DEFAULT_PR2_PATH)
                taxmap_path = None
            if err:
                accord_items.append(dbc.AccordionItem(html.P(err, style={"color": "red"}), title=f"{db_label} (error)"))
                db_errors.append(f"{db_label}: {err}")
                continue
            db_prefix = os.path.join(pr2_dir, "pr2")
            ok, log = ensure_blast_db(fasta_path, db_prefix)
            if not ok:
                accord_items.append(dbc.AccordionItem(html.Pre(log), title=f"{db_label} (makeblastdb failed)"))
                _write_blast_log(db_name, "", log, extra_notes=f"makeblastdb failed for {db_label}")
                db_errors.append(f"{db_label}: makeblastdb failed")
                continue
            stdout, stderr, code, attempts = _run_blastn_with_retry(query_fasta, db_prefix, output_path, max_hits=max_hits, evalue=evalue, remote=False)
            if code == 0 and taxmap_path:
                apply_pr2_taxmap_to_blast(output_path, taxmap_path)
        elif db_name == "custom":
            if not custom_db_path:
                accord_items.append(dbc.AccordionItem(html.P("Custom DB selected but no FASTA chosen.", style={"color": "red"}), title="CUSTOM (error)"))
                db_errors.append("CUSTOM: no FASTA chosen")
                continue
            db_prefix = os.path.join(os.path.dirname(custom_db_path), os.path.splitext(os.path.basename(custom_db_path))[0])
            ok, log = ensure_blast_db(custom_db_path, db_prefix)
            if not ok:
                accord_items.append(dbc.AccordionItem(html.Pre(log), title="CUSTOM (makeblastdb failed)"))
                _write_blast_log(db_name, "", log, extra_notes="makeblastdb failed for CUSTOM")
                db_errors.append("CUSTOM: makeblastdb failed")
                continue
            stdout, stderr, code, attempts = _run_blastn_with_retry(query_fasta, db_prefix, output_path, max_hits=max_hits, evalue=evalue, remote=False)
        else:
            continue

        if code != 0 or not os.path.exists(output_path):
            _write_blast_log(db_name, stdout, stderr, attempts=attempts if 'attempts' in locals() else None, extra_notes="BLAST command failed")
            db_errors.append(f"{db_label}: BLAST command failed")
            accord_items.append(dbc.AccordionItem([
                html.P("BLAST failed. Check logs.", style={"color": "red"}),
                html.H6("stdout"), html.Pre(stdout),
                html.H6("stderr"), html.Pre(stderr)
            ], title=f"{db_label} (failed)"))
            continue

        _write_blast_log(db_name, stdout, stderr, attempts=attempts if 'attempts' in locals() else None, extra_notes=f"BLAST {blast_version}")

        if _path_is_valid(output_path):
            db_outputs[f"blast_{db_name}"] = output_path
            _mirror_output_to_results(output_path, os.path.join(results_dir, f"blast_{db_name}.tsv"))

        preview = html.Div("No hits.")
        try:
            df = pd.read_csv(output_path, sep="\t", header=None, names=BLAST_TABLE_COLUMNS)
            preview = dash_table.DataTable(
                data=df.head(20).to_dict("records"),
                columns=[{"name": i, "id": i} for i in df.columns],
                page_size=10,
                style_table={"overflowX": "auto", "marginTop": "1em"},
                style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
                style_cell={"textAlign": "left", "padding": "5px"}
            )
        except Exception as e:
            preview = html.P(f"Could not parse BLAST output: {e}", style={"color": "orange"})

        accord_items.append(dbc.AccordionItem([
            html.P(db_method_info.get(db_name, ""), style={"color": "#6c757d"}) if db_name in db_method_info else html.Div(),
            html.A("Download BLAST results", href=_download_href(output_path), download=os.path.basename(output_path), target="_blank"),
            html.Hr(),
            preview
        ], title=f"{db_label} results"))

    comparison_df, pairwise_df = _build_blast_db_comparison(results_dir, db_filter=selected_dbs)
    if comparison_df is not None:
        comparison_path = os.path.join(db_outputs_dir, "blast_db_comparison.tsv")
        comparison_df.to_csv(comparison_path, sep="\t", index=False)
        if _path_is_valid(comparison_path):
            db_outputs["blast_db_comparison"] = comparison_path
            _mirror_output_to_results(comparison_path, os.path.join(results_dir, "blast_db_comparison.tsv"))

        comparison_preview = dash_table.DataTable(
            data=comparison_df.head(20).to_dict("records"),
            columns=[{"name": i, "id": i} for i in comparison_df.columns],
            page_size=10,
            style_table={"overflowX": "auto", "marginTop": "1em"},
            style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
            style_cell={"textAlign": "left", "padding": "5px"}
        )

        compare_children = [
            html.P("Comparison of top BLAST hit per query across selected databases using raw and normalized taxonomy labels."),
            html.A("Download cross-database comparison", href=_download_href(comparison_path), download=os.path.basename(comparison_path), target="_blank"),
            html.Hr(),
            comparison_preview
        ]

        if pairwise_df is not None and not pairwise_df.empty:
            pairwise_path = os.path.join(db_outputs_dir, "blast_db_pairwise_jaccard.tsv")
            pairwise_df.to_csv(pairwise_path, sep="\t", index=False)
            if _path_is_valid(pairwise_path):
                db_outputs["blast_db_pairwise_jaccard"] = pairwise_path
                _mirror_output_to_results(pairwise_path, os.path.join(results_dir, "blast_db_pairwise_jaccard.tsv"))
            pairwise_preview = dash_table.DataTable(
                data=pairwise_df.to_dict("records"),
                columns=[{"name": i, "id": i} for i in pairwise_df.columns],
                page_size=10,
                style_table={"overflowX": "auto", "marginTop": "1em"},
                style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
                style_cell={"textAlign": "left", "padding": "5px"}
            )
            compare_children.extend([
                html.Br(),
                html.P("Pairwise agreement is reported as exact label matching and normalized label matching for more robust cross-database interpretation."),
                html.A("Download pairwise DB agreement (Jaccard exact + normalized)", href=_download_href(pairwise_path), download=os.path.basename(pairwise_path), target="_blank"),
                html.Hr(),
                pairwise_preview
            ])

        accord_items.append(dbc.AccordionItem(compare_children, title="Cross-database comparison"))

    db_metrics_df = _build_db_metrics_table(results_dir, comparison_df=comparison_df, pairwise_df=pairwise_df, db_filter=selected_dbs)
    if db_metrics_df is not None and not db_metrics_df.empty:
        metrics_path = os.path.join(db_outputs_dir, "db_analysis_metrics.tsv")
        db_metrics_df.to_csv(metrics_path, sep="\t", index=False)
        if _path_is_valid(metrics_path):
            db_outputs["db_analysis_metrics"] = metrics_path
            _mirror_output_to_results(metrics_path, os.path.join(results_dir, "db_analysis_metrics.tsv"))
        metrics_preview = dash_table.DataTable(
            data=db_metrics_df.to_dict("records"),
            columns=[{"name": i, "id": i} for i in db_metrics_df.columns],
            page_size=10,
            style_table={"overflowX": "auto", "marginTop": "1em"},
            style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
            style_cell={"textAlign": "left", "padding": "5px"}
        )
        accord_items.append(dbc.AccordionItem([
            html.P("Consolidated metrics by selected database: hit coverage, taxonomic diversity, consensus support and pairwise agreement."),
            html.A("Download DB analysis metrics", href=_download_href(metrics_path), download=os.path.basename(metrics_path), target="_blank"),
            html.Hr(),
            metrics_preview
        ], title="Database metrics overview"))

    if not db_outputs:
        final_status = STEP_STATUS_SKIPPED if not accord_items else STEP_STATUS_FAILED
    else:
        final_status = STEP_STATUS_SUCCESS
    _write_db_manifest(
        final_status,
        outputs=db_outputs,
        error="; ".join(db_errors),
        notes="Database classification completed." if db_outputs else "No validated database outputs were produced.",
    )

    if not accord_items:
        return html.Div("No classification results were generated.")

    return dbc.Accordion(accord_items, start_collapsed=True, always_open=False, style={"marginTop": "1em"})
@app.callback(
    [Output("custom-db-status", "children"), Output("project-state", "data", allow_duplicate=True)],
    Input("custom-db-upload", "contents"),
    State("custom-db-upload", "filename"),
    State("project-state", "data"),
    prevent_initial_call=True
)
def upload_custom_db(contents, filename, project_state):
    if not contents or not filename:
        return dash.no_update, dash.no_update
    state = (project_state or {}).copy()
    pdir = state.get("project_dir")
    if not pdir:
        return html.Div("Create/select a project before uploading a custom DB.", style={"color": "orange"}), dash.no_update

    db_dir = os.path.join(pdir, "db", "custom")
    os.makedirs(db_dir, exist_ok=True)

    _, content_string = contents.split(",", 1)
    raw_path = os.path.join(db_dir, filename)
    with open(raw_path, "wb") as f:
        f.write(base64.b64decode(content_string))

    added = []
    if filename.lower().endswith(".zip"):
        extract_dir = os.path.join(db_dir, _safe_name(os.path.splitext(filename)[0]))
        os.makedirs(extract_dir, exist_ok=True)
        try:
            _safe_extract_zip(raw_path, extract_dir)
        except Exception as e:
            return html.Div(f"Custom DB ZIP extraction error: {e}", style={"color": "red"}), dash.no_update
        for root, _, files in os.walk(extract_dir):
            for f in files:
                if f.lower().endswith((".fasta", ".fa")):
                    added.append(os.path.join(root, f))
    else:
        if filename.lower().endswith((".fasta", ".fa")):
            added.append(raw_path)

    if not added:
        return html.Div("No FASTA files found in the custom DB upload.", style={"color": "orange"}), dash.no_update

    custom_dbs = state.get("custom_dbs", [])
    for p in added:
        if p not in custom_dbs:
            custom_dbs.append(p)
    state["custom_dbs"] = custom_dbs
    return html.Div(f"Custom DBs added: {len(added)}"), state

@app.callback(
    Output("custom-db-dropdown", "options"),
    Input("project-state", "data")
)
def refresh_custom_db_options(project_state):
    custom_dbs = (project_state or {}).get("custom_dbs", [])
    return [{"label": os.path.basename(p), "value": p} for p in custom_dbs]

@app.callback(
    Output('combined-fasta-status', 'children'),
    Input('build-combined-fasta-btn', 'n_clicks'),
    [State('samples-store', 'data'), State('project-state', 'data')],
    prevent_initial_call=True
)
def build_combined_fasta(n_clicks, samples_store, project_state):
    if not n_clicks:
        return dash.no_update
    samples = samples_store or {}
    if not isinstance(samples, dict) or not samples:
        return html.Div("No samples found. Upload or assign samples first.", style={"color": "orange"})
    fasta_list = []
    for _, info in samples.items():
        f = (info.get('outputs') or {}).get('input_fasta') or (info.get('paths') or {}).get('se_fasta')
        if f and os.path.exists(f) and os.path.getsize(f) > 0:
            fasta_list.append(f)
    if not fasta_list:
        return html.Div("No per-sample FASTA files available to combine.", style={"color": "orange"})
    results_dir = _resolve_output_path(project_state, "analysis", "final")
    out_path = os.path.join(results_dir, 'all_samples_combined.fasta')
    ok = write_combined_fasta(fasta_list, out_path)
    if not ok:
        return html.Div("Failed to build combined FASTA.", style={"color": "red"})
    link = html.A("Download combined FASTA", href=_download_href(out_path), target="_blank", download=os.path.basename(out_path))
    return html.Div([html.Span("Combined FASTA created: ", style={"marginRight": "8px"}), link])

@app.callback(
    Output("final-results", "children"),
    Input("run-completed-signal", "data"),
    State("project-state", "data"),
    State("run-config", "data"),
    State("workflow-scope-samples", "data"),
    prevent_initial_call=True
)
def build_final_results(run_completed_signal, project_state, run_config, workflow_scope):
    if not run_completed_signal:
        return ""

    results_dir = _resolve_output_path(project_state, "analysis", "final")
    project_dir = (project_state or {}).get("project_dir")
    run_dir = (project_state or {}).get("run_dir")

    analysis_manifest = _read_step_manifest(project_state, "analysis", default={}) or {}
    db_manifest = _read_step_manifest(project_state, "databases", default={}) or _read_step_manifest(project_state, "db", default={}) or {}
    phylo_manifest = _read_step_manifest(project_state, "phylo", default={}) or {}
    qc_manifest = _read_step_manifest(project_state, "qc", default={}) or {}
    preprocess_manifest = _read_step_manifest(project_state, "preprocess", default={}) or {}

    def _status_of(manifest, fallback_paths=None):
        status = (manifest or {}).get("status") if isinstance(manifest, dict) else None
        if status in {STEP_STATUS_SUCCESS, STEP_STATUS_FAILED, STEP_STATUS_SKIPPED, STEP_STATUS_RUNNING, STEP_STATUS_PENDING}:
            return status
        if fallback_paths and any(_path_is_valid(path) for path in fallback_paths):
            return STEP_STATUS_SUCCESS
        return STEP_STATUS_SKIPPED

    analysis_outputs = (analysis_manifest.get("outputs", {}) if isinstance(analysis_manifest, dict) else {}) or {}
    analysis_artifacts = (analysis_outputs.get("artifacts", {}) if isinstance(analysis_outputs, dict) else {}) or {}
    db_outputs = (db_manifest.get("outputs", {}) if isinstance(db_manifest, dict) else {}) or {}
    phylo_outputs = (phylo_manifest.get("outputs", {}) if isinstance(phylo_manifest, dict) else {}) or {}

    otu_representatives_path = _manifest_path_value(analysis_outputs, "vsearch", "outputs", "otus_fasta") or _manifest_path_value(analysis_outputs, "vsearch", "otus_fasta") or os.path.join(results_dir, "otus.fasta")
    asv_fasta_path = _manifest_path_value(analysis_outputs, "dada2", "outputs", "asv_fasta") or _manifest_path_value(analysis_outputs, "dada2", "asv_fasta") or os.path.join(results_dir, "asv.fasta")
    otu_abundance_matrix_path = _manifest_path_value(analysis_artifacts, "otu_sample_matrix") or _manifest_path_value(analysis_artifacts, "asv_sample_matrix") or os.path.join(results_dir, "otu_sample_matrix.tsv")
    asv_abundance_matrix_path = _manifest_path_value(analysis_artifacts, "asv_sample_matrix") or os.path.join(results_dir, "asv_sample_matrix.tsv")
    alpha_diversity_path = _manifest_path_value(analysis_artifacts, "alpha_diversity") or os.path.join(results_dir, "alpha_diversity.tsv")
    beta_diversity_dir = os.path.join(results_dir, "beta_diversity")
    comparative_dir = os.path.join(results_dir, "comparative_analysis")
    qc_summary_path = os.path.join(results_dir, "qc_summary.tsv")
    summary_stats_path = _manifest_path_value(analysis_artifacts, "summary_stats") or os.path.join(results_dir, "summary_stats.tsv")
    stacked_abundance_path = _manifest_path_value(analysis_artifacts, "stacked_abundance") or os.path.join(comparative_dir, "stacked_abundance.tsv")
    stacked_abundance_html_path = _manifest_path_value(analysis_artifacts, "stacked_abundance_html") or os.path.join(comparative_dir, "stacked_abundance.html")
    rarefaction_curves_path = _manifest_path_value(analysis_artifacts, "rarefaction_curves") or os.path.join(comparative_dir, "rarefaction_curves.tsv")
    rarefaction_curves_html_path = _manifest_path_value(analysis_artifacts, "rarefaction_curves_html") or os.path.join(comparative_dir, "rarefaction_curves.html")
    blast_summary_path = os.path.join(results_dir, "blast_summary.tsv")
    blast_compare_path = os.path.join(results_dir, "blast_db_comparison.tsv")
    blast_pairwise_path = os.path.join(results_dir, "blast_db_pairwise_jaccard.tsv")
    db_metrics_path = os.path.join(results_dir, "db_analysis_metrics.tsv")
    analysis_overview_path = os.path.join(results_dir, "analysis_overview_comparison.tsv")
    jaccard_distance_path = os.path.join(beta_diversity_dir, "jaccard_distance.tsv")
    braycurtis_distance_path = os.path.join(beta_diversity_dir, "braycurtis_distance.tsv")
    pcoa_coords_path = os.path.join(beta_diversity_dir, "pcoa_coords.tsv")
    nmds_coords_path = os.path.join(beta_diversity_dir, "nmds_coords.tsv")
    shared_taxa_path = os.path.join(comparative_dir, "shared_taxa.tsv")
    unique_taxa_path = os.path.join(comparative_dir, "unique_taxa_per_sample.tsv")
    heatmap_abundance_path = os.path.join(comparative_dir, "heatmap_abundance.tsv")
    mafft_alignment_path = _manifest_path_value(phylo_outputs, "aligned_fasta") or _manifest_path_value(phylo_outputs, "aligned_safe_ids") or os.path.join(results_dir, "phylo", "aligned_safe_ids.fasta")
    fasttree_tree_path = _manifest_path_value(phylo_outputs, "fasttree_tree") or _manifest_path_value(phylo_outputs, "fasttree_newick") or os.path.join(results_dir, "phylo", "tree.nwk")
    raxml_best_tree_path = _manifest_path_value(phylo_outputs, "raxml_best_tree") or os.path.join(_resolve_output_path(project_state, "phylo", "final"), "raxml.raxml.bestTree")

    step_statuses = {
        "qc": _status_of(qc_manifest, [qc_summary_path]),
        "preprocess": _status_of(preprocess_manifest),
        "analysis": _status_of(analysis_manifest, [otu_representatives_path, asv_fasta_path, otu_abundance_matrix_path]),
        "databases": _status_of(db_manifest, [blast_summary_path, blast_compare_path, blast_pairwise_path, db_metrics_path]),
        "phylo": _status_of(phylo_manifest, [fasttree_tree_path, raxml_best_tree_path]),
    }

    def _add_link(links, label, path):
        if _path_is_valid(path):
            links.append(html.Li(html.A(label, href=_download_href(path), target="_blank", download=os.path.basename(path))))

    links = []
    _add_link(links, "QC summary", qc_summary_path)
    _add_link(links, "OTU representatives FASTA", otu_representatives_path)
    _add_link(links, "ASV FASTA", asv_fasta_path)
    _add_link(links, "OTU abundance matrix", otu_abundance_matrix_path)
    _add_link(links, "ASV abundance matrix", asv_abundance_matrix_path)
    _add_link(links, "Alpha diversity", alpha_diversity_path)
    _add_link(links, "Summary stats", summary_stats_path)
    _add_link(links, "Jaccard beta diversity", jaccard_distance_path)
    _add_link(links, "Bray-Curtis beta diversity", braycurtis_distance_path)
    _add_link(links, "PCoA coordinates", pcoa_coords_path)
    _add_link(links, "NMDS coordinates", nmds_coords_path)
    _add_link(links, "Shared taxa", shared_taxa_path)
    _add_link(links, "Unique taxa per sample", unique_taxa_path)
    _add_link(links, "Heatmap abundance", heatmap_abundance_path)
    _add_link(links, "Stacked abundance", stacked_abundance_path)
    _add_link(links, "Stacked abundance HTML", stacked_abundance_html_path)
    _add_link(links, "Rarefaction curves", rarefaction_curves_path)
    _add_link(links, "Rarefaction HTML", rarefaction_curves_html_path)
    _add_link(links, "BLAST summary", blast_summary_path)
    _add_link(links, "BLAST cross-database comparison", blast_compare_path)
    _add_link(links, "BLAST pairwise DB agreement", blast_pairwise_path)
    _add_link(links, "DB analysis metrics", db_metrics_path)
    _add_link(links, "Analysis overview comparison", analysis_overview_path)
    _add_link(links, "MAFFT alignment", mafft_alignment_path)
    _add_link(links, "FastTree Newick", fasttree_tree_path)
    _add_link(links, "RAxML best tree", raxml_best_tree_path)
    _add_link(links, "Analysis manifest", _resolve_output_path(project_state, "analysis", "final", "analysis_manifest.json"))
    _add_link(links, "Final manifest", _final_manifest_path(project_state))

    summary_note = []
    for step_name, status in step_statuses.items():
        summary_note.append(f"{step_name}={status}")

    blast_summary = None
    compare_df = None
    pairwise_df = None
    db_metrics_df = None
    analysis_overview_df = None
    db_classification_dir = _resolve_output_path(project_state, "db", "final", "db_classification")
    if _status_of(db_manifest, [blast_summary_path, blast_compare_path, blast_pairwise_path, db_metrics_path]) == STEP_STATUS_SUCCESS:
        if _path_is_valid(blast_summary_path):
            try:
                blast_summary = pd.read_csv(blast_summary_path, sep="\t", header=None, names=BLAST_TABLE_COLUMNS)
            except Exception:
                blast_summary = None
        if _path_is_valid(blast_compare_path):
            try:
                compare_df = pd.read_csv(blast_compare_path, sep="\t")
            except Exception:
                compare_df = None
        if _path_is_valid(blast_pairwise_path):
            try:
                pairwise_df = pd.read_csv(blast_pairwise_path, sep="\t")
            except Exception:
                pairwise_df = None
        if _path_is_valid(db_metrics_path):
            try:
                db_metrics_df = pd.read_csv(db_metrics_path, sep="\t")
            except Exception:
                db_metrics_df = None
        if _path_is_valid(analysis_overview_path):
            try:
                analysis_overview_df = pd.read_csv(analysis_overview_path, sep="\t")
            except Exception:
                analysis_overview_df = None

    def _table_preview(df, table_id):
        if df is None or df.empty:
            return html.Div("No rows available for preview.", style={"color": "#6c757d"})
        return dash_table.DataTable(
            id=table_id,
            data=df.head(20).to_dict("records"),
            columns=[{"name": i, "id": i} for i in df.columns],
            page_size=10,
            style_table={"overflowX": "auto", "marginTop": "0.5em"},
            style_header={"backgroundColor": "rgb(230, 230, 230)", "fontWeight": "bold"},
            style_cell={"textAlign": "left", "padding": "5px"}
        )

    db_previews = []
    if _status_of(db_manifest, [blast_summary_path, blast_compare_path, blast_pairwise_path, db_metrics_path]) == STEP_STATUS_SUCCESS:
        for i, blast_file in enumerate(_select_blast_files(db_classification_dir), start=1):
            try:
                df_blast = pd.read_csv(blast_file, sep="\t", header=None, names=BLAST_TABLE_COLUMNS)
                db_previews.append(html.Div([
                    html.H6(f"{os.path.basename(blast_file)} (top 20 rows)"),
                    _table_preview(df_blast, f"final-db-blast-preview-{i}")
                ]))
            except Exception as e:
                db_previews.append(html.Div(f"Could not preview {os.path.basename(blast_file)}: {e}", style={"color": "orange"}))

    summary_previews = []
    if blast_summary is not None:
        summary_previews.append(html.Div([html.H6("BLAST summary (top 20 rows)"), _table_preview(blast_summary, "final-blast-summary-preview")]))
    if compare_df is not None:
        summary_previews.append(html.Div([html.H6("Cross-database comparison (top 20 rows)"), _table_preview(compare_df, "final-blast-compare-preview")]))
    if pairwise_df is not None and not pairwise_df.empty:
        summary_previews.append(html.Div([html.H6("Pairwise database agreement (top 20 rows)"), _table_preview(pairwise_df, "final-blast-pairwise-preview")]))
    if db_metrics_df is not None and not db_metrics_df.empty:
        summary_previews.append(html.Div([html.H6("Database metrics overview (top 20 rows)"), _table_preview(db_metrics_df, "final-db-metrics-preview")]))
    if analysis_overview_df is not None and not analysis_overview_df.empty:
        summary_previews.append(html.Div([html.H6("Analysis overview comparison (top 20 rows)"), _table_preview(analysis_overview_df, "final-analysis-overview-preview")]))

    final_manifest = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project_name": (project_state or {}).get("project_name", ""),
        "run_dir": run_dir or "",
        "step_statuses": step_statuses,
        "outputs": {path_label: path for path_label, path in {
            "qc_summary": qc_summary_path,
            "otu_representatives": otu_representatives_path,
            "asv_fasta": asv_fasta_path,
            "otu_abundance_matrix": otu_abundance_matrix_path,
            "asv_abundance_matrix": asv_abundance_matrix_path,
            "alpha_diversity": alpha_diversity_path,
            "summary_stats": summary_stats_path,
            "jaccard_distance": jaccard_distance_path,
            "braycurtis_distance": braycurtis_distance_path,
            "pcoa_coords": pcoa_coords_path,
            "nmds_coords": nmds_coords_path,
            "shared_taxa": shared_taxa_path,
            "unique_taxa_per_sample": unique_taxa_path,
            "heatmap_abundance": heatmap_abundance_path,
            "heatmap_abundance_html": _manifest_path_value(analysis_artifacts, "heatmap_abundance_html") or os.path.join(comparative_dir, "heatmap_abundance.html"),
            "stacked_abundance": stacked_abundance_path,
            "stacked_abundance_html": stacked_abundance_html_path,
            "rarefaction_curves": rarefaction_curves_path,
            "rarefaction_curves_html": rarefaction_curves_html_path,
            "blast_summary": blast_summary_path,
            "blast_db_comparison": blast_compare_path,
            "blast_db_pairwise_jaccard": blast_pairwise_path,
            "db_analysis_metrics": db_metrics_path,
            "analysis_overview_comparison": analysis_overview_path,
            "mafft_alignment": mafft_alignment_path,
            "fasttree_newick": fasttree_tree_path,
            "raxml_best_tree": raxml_best_tree_path,
        }.items() if _path_is_valid(path)},
        "workflow_scope": workflow_scope or {},
    }
    _write_json(_final_manifest_path(project_state), final_manifest)

    run_size_bytes = _get_directory_size(run_dir) if run_dir else 0
    if run_dir:
        session_json_path = _run_session_path(run_dir)
        session_payload = _read_json(session_json_path, default={}) or {}
        if not isinstance(session_payload, dict):
            session_payload = {}
        session_payload["completed_at"] = datetime.now().isoformat(timespec="seconds")
        session_payload["storage_bytes"] = run_size_bytes
        session_payload.setdefault("outputs", {})
        session_payload["outputs"].update(final_manifest["outputs"])
        _write_json(session_json_path, session_payload)
        if DELETE_AFTER_STEPS and not KEEP_INTERMEDIATES:
            _apply_cleanup_policy(run_dir)

    session_zip = _create_session_results_zip(project_state)
    if project_dir and run_dir:
        _apply_cleanup_post_run(project_dir, run_dir)

    if session_zip and _path_is_valid(session_zip):
        zip_name = os.path.basename(session_zip)
        zip_size = os.path.getsize(session_zip) / (1024 * 1024)
        links.append(html.Li(html.A(f"Download {zip_name} ({zip_size:.1f} MB)", href=_download_href(session_zip), target="_blank", download=zip_name)))

    export_note = html.Div([
        html.Hr(),
        html.H6("Run Export"),
        html.Ul([
            html.Li(
                html.A(
                    f"Download {os.path.basename(session_zip)} ({os.path.getsize(session_zip) / (1024 * 1024):.1f} MB)",
                    href=_download_href(session_zip),
                    target="_blank",
                    download=os.path.basename(session_zip),
                )
            )
        ]) if session_zip and _path_is_valid(session_zip) else html.Div(
            _LAST_EXPORT_ERROR or "No export ZIP was generated.",
            style={"color": "orange"}
        )
    ])

    selected_steps = (run_config or {}).get("steps", []) if isinstance(run_config, dict) else []
    visible_statuses = {}
    for step, status in step_statuses.items():
        if step == "databases" and "db" in selected_steps:
            visible_statuses[step] = status
        elif step in selected_steps:
            visible_statuses[step] = status
    
    result_children = [
        html.H5("Final Results"),
        html.P("; ".join(summary_note)),
        html.H6("Step Statuses"),
        html.Ul([html.Li(f"{step}: {status}") for step, status in visible_statuses.items()]) if visible_statuses else html.Div("No selected step statuses to display.", style={"color": "#6c757d"}),
        html.Ul(links) if links else html.Div("No validated outputs were produced.", style={"color": "orange"}),
    ]

    if "db" in selected_steps:
        result_children.append(html.Hr())
        result_children.append(html.H5("Database Outputs Preview"))
        result_children.append(html.Div(db_previews if db_previews else [html.Div("No database step executed or no validated BLAST tables found.", style={"color": "#6c757d"})]))
        result_children.append(html.Hr())
        result_children.append(html.H5("Database Summary Tables Preview"))
        result_children.append(html.Div(summary_previews if summary_previews else [html.Div("No consolidated database summary tables found for preview.", style={"color": "#6c757d"})]))

    result_children.append(export_note)
    return html.Div(result_children)

@app.callback(
    Output("theme-root", "className"),
    Input("theme-toggle", "value"),
    prevent_initial_call=False
)
def set_theme(is_dark):
    mode = "dark" if is_dark else "light"
    return f"theme-{mode}"

if __name__ == '__main__':
    app.run(debug=True)