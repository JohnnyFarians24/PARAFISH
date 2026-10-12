"""Orchestration callbacks for running the full pipeline from the Setup page."""

from __future__ import annotations

import json
import threading
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO
from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from app import _pipeline_threads, _progress_lock
from analysis.analysis_metrics import build_metrics_from_outputs, feature_sets_from_outputs
from analysis.dada2_runner import run_dada2
from analysis.feature_set_comparison import build_feature_set_comparison
from analysis.vsearch_runner import run_vsearch_filter_sample, run_vsearch_pool_and_cluster
from core.config import MIN_UNIQUE_SIZE, NCBI_USE_LOCAL_DB, get_db_paths
from core.preflight import preflight_notes
from core.progress_writer import append_log, read_progress, write_progress
from core.storage import _apply_cleanup_post_run
from databases.bold_db import BOLDConfig
from databases.database_metrics import build_database_metrics
from databases.ncbi_db import NCBIConfig
from databases.pr2_db import PR2Config
from databases.reference_filter import marker_filter_applies, resolve_run_marker
from databases.silva_db import SILVAConfig
from databases.taxonomy_runner import run_taxonomy_pipeline
from phylo import phylo_input_error
from phylo.phylo_metrics import build_metrics_from_files
from phylo.phylo_runner import run_phylo as run_phylo_step
from preprocess.cutadapt_runner import resolve_primers, run_cutadapt
from preprocess.pear_runner import run_pear
from preprocess.preprocess_metrics import build_preprocess_metrics
from preprocess.trimmomatic_runner import run_trimmomatic
from qc.fastqc_runner import run_fastqc
from qc.qc_metrics import build_qc_metrics
from ui.pipeline_state import stamp_record, utc_now_iso
from ui.runtime_state import build_sample_from_record, run_config_from_store, serialize_step_output

def _stage_entry(stage_name: str, status: str) -> str:
    return f"{stage_name}: {status}"


def _log_entry(ts: str, level: str, msg: str) -> dict:
    return {"ts": ts, "level": level, "msg": msg}


def _failed_record(run_id: str, outputs: dict, logs: list) -> dict:
    return {
        "signature": run_id,
        "status": "failure",
        "updated_at": utc_now_iso(),
        "outputs": outputs,
        "log_lines": logs.copy(),
    }


def _count_features(fasta_path) -> int:
    try:
        with Path(fasta_path).open("r", errors="replace") as handle:
            return sum(1 for line in handle if line.startswith(">"))
    except OSError:
        return 0


def _selected_samples(samples):
	return [sample for sample in samples or [] if getattr(sample, "selected", True)]


def _pipeline_run_marker(run_config) -> str:
    """The marker this run amplified, resolved the way Setup resolves it.

    Primer preset first, the samples' file-name markers as the fallback — the same
    ``resolve_run_marker`` the Setup page, the standalone Databases callback and
    ``RunConfig.marker`` already use. This thread used ``resolve_marker`` over the file
    names alone, so a COI/Leray run on files named ``aptl64b_01`` reached the Databases
    stage as OTHER and SILVA/PR2 were not marked not_applicable. With no preset (or
    ``auto``) the preset names no marker and the result is the file-name marker, as before.
    """
    marker, _ = resolve_run_marker((sample.marker for sample in run_config.samples), run_config.vsearch_params)
    return marker


def _recorded_parameters(run_config, requested_sources, ncbi_local: bool) -> dict:
    """Every stage's parameters as this run will use them, for the progress.json record.

    The run_config block used to name the sources and little else, so the clustering
    identity, the SINTAX cutoff or the bootstrap count of a finished run could only be
    recovered from tool logs. The ``databases`` values are computed with the SAME
    expressions the Databases stage uses further down in ``_run_pipeline_thread`` (the
    NCBIConfig/BOLDConfig construction) — if one changes, change the other.
    """
    from databases.taxonomy_runner import SINTAX_RANDSEED

    vsearch_params = dict(run_config.vsearch_params or {})
    databases = {
        "db_sources": list(requested_sources),
        "marker": _pipeline_run_marker(run_config),
        "sintax_cutoff": float(vsearch_params.get("sintax_cutoff", 0.8)),
        "sintax_randseed": SINTAX_RANDSEED,
        "consensus_method": str(vsearch_params.get("db_consensus") or "lca"),
    }
    if "ncbi" in requested_sources:
        ncbi_defaults = NCBIConfig.__dataclass_fields__
        databases["ncbi"] = {
            "mode": "local" if ncbi_local else "remote",
            "task": str(vsearch_params.get("ncbi_task") or ("blastn" if ncbi_local else "megablast")),
            "word_size": int(vsearch_params["ncbi_word_size"]) if vsearch_params.get("ncbi_word_size") else (11 if ncbi_local else None),
            "evalue": float(vsearch_params.get("ncbi_evalue", 1e-10)),
            "max_target_seqs": int(vsearch_params.get("ncbi_maxhits", 20)),
            "perc_identity": ncbi_defaults["perc_identity"].default,
            "qcov_hsp_perc": ncbi_defaults["qcov_hsp_perc"].default,
        }
    if "bold" in requested_sources:
        bold_params = run_config.bold_params or {}
        databases["bold"] = {
            "db": int(bold_params.get("db", 1)),
            "mode": int(bold_params.get("mode", 2)),
            "thresholds": list(BOLDConfig.__dataclass_fields__["thresholds"].default_factory()),
        }
    phylo_mode = getattr(run_config, "phylo_mode", None)
    return {
        "enabled_stages": dict(run_config.enabled_stages or {}),
        "trimmomatic": dict(run_config.trimmomatic_params or {}),
        "pear": dict(run_config.pear_params or {}),
        "vsearch": vsearch_params,
        "dada2": dict(run_config.dada2_params or {}),
        "databases": databases,
        "phylo": {
            "mode": getattr(phylo_mode, "value", phylo_mode),
            "mafft": dict(run_config.mafft_params or {}),
            "tree": dict(run_config.phylo_params or {}),
        },
    }


def _probe_version(argv, pattern: str) -> str:
    """One tool's version from its own output, or ``"unresolved"``. Never raises."""
    import re
    import subprocess

    try:
        out = subprocess.run(argv, capture_output=True, text=True, timeout=15, input="")
        match = re.search(pattern, f"{out.stdout or ''}\n{out.stderr or ''}")
        return match.group(1) if match else "unresolved"
    except Exception:
        return "unresolved"


def _recorded_tool_versions(run_config):
    """Tool -> version, resolved NOW from the configured binaries. Never raises.

    PROVENANCE.txt used to resolve versions when the ZIP was exported, so an export made
    after an upgrade reported a version that never touched the data. The resolver is the
    one the Results page already uses (configured paths, not $PATH), plus the tools it
    does not list. A tool that cannot be resolved is ``"unresolved"`` — the run goes on.
    """
    try:
        from ui.callbacks_results_summary import _resolve_tool_versions

        versions = {}
        for tool, version, resolved in _resolve_tool_versions():
            detected = bool(resolved) and "not detected" not in str(version)
            versions[str(tool)] = str(version) if detected else "unresolved"
    except Exception:
        return "unresolved"

    try:
        import shutil

        from core.config import CUTADAPT_PATH, RSCRIPT_PATH, TRIMMOMATIC_JAR_PATH

        versions["cutadapt"] = _probe_version([CUTADAPT_PATH, "--version"], r"^\s*([\d.]+)") if shutil.which(CUTADAPT_PATH) or Path(CUTADAPT_PATH).is_file() else "unresolved"
        versions["Trimmomatic"] = (
            _probe_version(["java", "-jar", str(TRIMMOMATIC_JAR_PATH), "-version"], r"([\d.]+)")
            if Path(TRIMMOMATIC_JAR_PATH).is_file() and shutil.which("java")
            else "unresolved"
        )
        if getattr(run_config, "analysis_method", "vsearch") in ("dada2", "both"):
            versions["DADA2"] = (
                _probe_version([RSCRIPT_PATH, "-e", "cat(as.character(packageVersion('dada2')))"], r"([\d.]+)")
                if shutil.which(RSCRIPT_PATH) or Path(RSCRIPT_PATH).is_file()
                else "unresolved"
            )
        try:
            from importlib.metadata import version as _package_version

            versions["boldigger3"] = _package_version("boldigger3")
        except Exception:
            versions["boldigger3"] = "unresolved"
    except Exception:
        pass
    return versions


def _recorded_references(run_config, requested_sources, db_paths, ncbi_local: bool) -> list:
    """Name, size, mtime and sha256 of each LOCAL reference this run will use.

    Only when the Databases stage is enabled, and only for the sources requested: a run
    that never opens SILVA has no business hashing 800 MB of it. Remote sources (NCBI nt,
    BOLD) have no file to fingerprint; their settings are under ``parameters.databases``.
    """
    from core.export import reference_fingerprint

    if not (run_config.enabled_stages or {}).get("databases", True):
        return []
    candidates = []
    if "silva" in requested_sources:
        candidates.append(("silva", getattr(db_paths, "silva", None)))
    if "pr2" in requested_sources:
        candidates.append(("pr2", getattr(db_paths, "pr2", None)))
        candidates.append(("pr2_taxonomy", getattr(db_paths, "pr2_tax", None)))
    if "ncbi" in requested_sources and ncbi_local:
        candidates.append(("ncbi_local", getattr(db_paths, "blast", None)))
    custom = (run_config.vsearch_params or {}).get("custom_db_path")
    if custom:
        candidates.append(("custom", custom))
    references = []
    for source, path in candidates:
        if not path:
            continue
        # A missing file is still listed (sha256 "unavailable"), except for the references
        # the Databases stage itself drops silently when the file is absent.
        if source != "ncbi_local" and not Path(path).exists():
            continue
        references.append({"source": source, **reference_fingerprint(path)})
    return references


def _run_provenance(run_config, run_config_data) -> dict:
    """The resolved config a run was made under, for the progress.json record.

    progress.json used to carry only per-stage status, so two runs that differ only
    in configuration — BOLD on vs off, a different min_unique_size, primers stripped
    or not — produced different outputs that read as a regression against golden when
    they are nothing of the sort. This block makes the difference legible: a golden
    diff can be read under a known config instead of guessing.

    Everything here is resolvable before any stage runs (it is written at init, so it
    survives even a run whose tools fail), which is why it records the *requested*
    db_sources rather than the post-existence-filter set the Databases stage computes.
    """
    data = run_config_data or {}
    vsearch_params = run_config.vsearch_params or {}
    requested_sources = [
        str(source).lower()
        for source in (data.get("db_sources") or vsearch_params.get("db_sources") or [])
    ]
    forward, reverse = resolve_primers(run_config)
    markers = sorted({sample.marker.value for sample in _selected_samples(run_config.samples)})
    db_paths = get_db_paths()
    ncbi_local = bool(NCBI_USE_LOCAL_DB and db_paths.blast and db_paths.blast.exists())
    block = {
        "db_sources": requested_sources,
        "min_unique_size": int(vsearch_params.get("min_unique_size", MIN_UNIQUE_SIZE)),
        "analysis_method": getattr(run_config, "analysis_method", "vsearch"),
        "branch": run_config.branch.value,
        "markers": markers,
        "primers": {"fwd": forward, "rev": reverse},
        "ncbi_mode": "local" if ncbi_local else "remote",
    }
    # Recorded when the run starts, and ADDED to the keys above — never instead of them.
    # Each is best-effort on its own: describing the run must not be able to fail it, and
    # one that cannot be completed says so rather than taking the others down with it.
    # ``default=str`` because a parameter dict may carry a Path, which json cannot write.
    try:
        block["parameters"] = json.loads(json.dumps(_recorded_parameters(run_config, requested_sources, ncbi_local), default=str))
    except Exception:
        block["parameters"] = "unavailable"
    try:
        block["tool_versions"] = _recorded_tool_versions(run_config)
    except Exception:
        block["tool_versions"] = "unresolved"
    try:
        block["references"] = _recorded_references(run_config, requested_sources, db_paths, ncbi_local)
    except Exception:
        block["references"] = "unavailable"
    try:
        from core.export import git_commit

        block["parafish_commit"] = git_commit()
    except Exception:
        block["parafish_commit"] = "unknown"
    return block


def _run_pipeline_thread(run_config, run_config_data, progress_path):
    """Runs in a background thread and writes progress to disk."""
    stages: dict[str, str] = {}
    logs: list[dict] = []
    active_stage_results: list[bool] = []
    progress_state = {
        "run_id": run_config.run_id,
        "run_config": _run_provenance(run_config, run_config_data),
        "qc_record": {"status": "pending"},
        "preprocess_record": {"status": "pending"},
        "analysis_record": {"status": "pending"},
        "db_record": {"status": "pending"},
        "phylo_record": {"status": "pending"},
        "log_lines": [],
    }
    progress_path.parent.mkdir(parents=True, exist_ok=True)
    progress_path.write_text(json.dumps(progress_state), encoding="utf-8")

    qc_record = {"signature": run_config.run_id, "status": "skipped", "updated_at": utc_now_iso(), "outputs": {}}
    preprocess_record = {"signature": run_config.run_id, "status": "skipped", "updated_at": utc_now_iso(), "outputs": {}, "next_samples": []}
    analysis_record = {"signature": run_config.run_id, "status": "skipped", "updated_at": utc_now_iso(), "outputs": {}}
    db_record = {"signature": run_config.run_id, "status": "skipped", "updated_at": utc_now_iso(), "outputs": {}}
    phylo_record = {"signature": run_config.run_id, "status": "skipped", "updated_at": utc_now_iso(), "outputs": {}}
    # Bound here, not in the Preprocess branch: the Analysis stage refreshes the
    # read-tracking waterfall, and Analysis can run on a Preprocess that was skipped
    # (re-analysing an existing run).
    preprocess_metrics_files: dict[str, str] = {}

    def _publish() -> None:
        # Persist the stage badges the instant they change. write_progress merges the
        # "stages" dict (core/progress_writer), so this never clobbers a record,
        # run_config or log_lines — it only lets the poller see a stage as "running"
        # the moment it starts, instead of one stage behind.
        write_progress(progress_path, {"stages": dict(stages)})

    def _write_progress_state() -> None:
        # Note: "log_lines" is deliberately NOT included here — append_log() already
        # persists each entry incrementally, and write_progress()'s merge semantics
        # for "log_lines" concatenate, so re-sending the full in-memory history on
        # every stage checkpoint would duplicate everything already on disk.
        progress_state.update({
            "stages": dict(stages),
        })
        write_progress(progress_path, progress_state)

    def _skipped_record() -> dict:
        return {"status": "skipped", "outputs": {}, "updated_at": utc_now_iso()}

    def _ensure_existing_path(path_value, stage_name: str, sample_id: str, label: str):
        if path_value is None or path_value == "":
            entry = _log_entry(utc_now_iso(), "ERR", f"{stage_name}: {label} missing for {sample_id}")
            logs.append(entry)
            append_log(progress_path, entry)
            return None
        candidate = Path(path_value)
        if not candidate.exists():
            entry = _log_entry(utc_now_iso(), "ERR", f"{stage_name}: {label} does not exist for {sample_id}: {candidate}")
            logs.append(entry)
            append_log(progress_path, entry)
            return None
        return candidate

    try:
        enabled = run_config.enabled_stages or {}
        run_qc = enabled.get("qc", True)
        run_preprocess = enabled.get("preprocess", True)
        run_analysis = enabled.get("analysis", True)
        run_databases = enabled.get("databases", True)
        run_phylo_stage = enabled.get("phylo", True)
        analysis_method = getattr(run_config, "analysis_method", "vsearch")
        run_vsearch_branch = analysis_method in ("vsearch", "both")
        run_dada2_branch = analysis_method in ("dada2", "both")

        # QC
        qc_outputs = {}
        if run_qc:
            stages["QC"] = "running"
            _publish()
            entry = _log_entry(utc_now_iso(), "RUN", "Starting QC stage")
            logs.append(entry)
            append_log(progress_path, entry)
            all_ok = True
            for sample in run_config.samples:
                if _ensure_existing_path(sample.paths.r1, "QC", sample.internal_id, "r1") is None:
                    all_ok = False
                    continue
                if sample.mode.value == "PE" and _ensure_existing_path(sample.paths.r2, "QC", sample.internal_id, "r2") is None:
                    all_ok = False
                    continue
                res = run_fastqc(sample, run_config.output_dir / "qc" / sample.internal_id, threads=int(run_config.vsearch_params.get("threads", 4)))
                qc_outputs[sample.internal_id] = serialize_step_output(res)
                entry = _log_entry(
                    utc_now_iso(),
                    "OK" if res.success else "ERR",
                    f"FastQC {sample.internal_id}: {res.error_message or 'OK'}" if not res.success else f"FastQC {sample.internal_id}: OK",
                )
                logs.append(entry)
                append_log(progress_path, entry)
                if not res.success:
                    all_ok = False
            # Compute the QC tables once, here, so the UI callbacks only read CSVs.
            # Kept beside "outputs" rather than inside it: "outputs" is keyed by
            # sample id, and an extra key there would be read back as a sample.
            qc_metrics_files: dict[str, str] = {}
            try:
                qc_metrics_files = {
                    key: str(path) for key, path in build_qc_metrics(run_config.output_dir / "qc").items()
                }
            except Exception as exc:
                entry = _log_entry(utc_now_iso(), "ERR", f"Failed to build QC metrics: {exc}")
                logs.append(entry)
                append_log(progress_path, entry)
            qc_record = {
                "signature": run_config.run_id,
                "status": "success" if all_ok else "failure",
                "updated_at": utc_now_iso(),
                "outputs": qc_outputs,
                "metrics_files": qc_metrics_files,
                "log_lines": logs.copy(),
            }
            progress_state["qc_record"] = qc_record
            # Mark the stage terminal BEFORE the write, so the record and its badge
            # land in the same snapshot — never "success record, running badge".
            stages["QC"] = qc_record["status"]
            _write_progress_state()
            active_stage_results.append(qc_record["status"] == "success")
            _publish()
        else:
            qc_record = _skipped_record()
            progress_state["qc_record"] = qc_record
            stages["QC"] = "skipped"
            _write_progress_state()
            _publish()

        # Preprocess
        preprocess_outputs = {}
        next_samples = []
        if run_preprocess:
            stages["Preprocess"] = "running"
            _publish()
            entry = _log_entry(utc_now_iso(), "RUN", "Starting Preprocess stage")
            logs.append(entry)
            append_log(progress_path, entry)
            all_ok = True
            for sample in run_config.samples:
                entry = _log_entry(utc_now_iso(), "DBG", f"Preprocess sample: id={sample.internal_id}, r1={sample.paths.r1}, exists={Path(sample.paths.r1).exists() if sample.paths.r1 else False}")
                logs.append(entry)
                append_log(progress_path, entry)
                if _ensure_existing_path(sample.paths.r1, "Preprocess", sample.internal_id, "r1") is None:
                    all_ok = False
                    continue
                if sample.mode.value == "PE" and _ensure_existing_path(sample.paths.r2, "Preprocess", sample.internal_id, "r2") is None:
                    all_ok = False
                    continue
                # Remove PCR primers on the raw reads before quality trimming. When no
                # primers are configured (or cutadapt is missing) this passes the reads
                # through untouched, so behaviour is unchanged for runs without primers.
                from core.models import SamplePaths

                cut_res = run_cutadapt(sample, run_config)
                preprocess_outputs.setdefault(sample.internal_id, {})["cutadapt"] = serialize_step_output(cut_res)
                if cut_res.error_message:
                    entry = _log_entry(utc_now_iso(), "WARN", f"cutadapt {sample.internal_id}: {cut_res.error_message}")
                    logs.append(entry)
                    append_log(progress_path, entry)
                sample_for_trim = sample
                if (cut_res.metrics or {}).get("primers_removed"):
                    # Sequence-aware removal happened here, so the Analysis stage must not
                    # strip again by length (see analysis.vsearch_runner._primer_strip_lengths).
                    run_config.vsearch_params["primers_removed_upstream"] = True
                    r1_primer = cut_res.files.get("r1_primer")
                    r2_primer = cut_res.files.get("r2_primer")
                    sample_for_trim = replace(
                        sample,
                        paths=SamplePaths(r1=Path(r1_primer), r2=Path(r2_primer) if r2_primer else None),
                    )
                    entry = _log_entry(utc_now_iso(), "OK", f"cutadapt {sample.internal_id}: primers removidos")
                    logs.append(entry)
                    append_log(progress_path, entry)
                trim_res = run_trimmomatic(sample_for_trim, run_config)
                preprocess_outputs.setdefault(sample.internal_id, {})["trimmomatic"] = serialize_step_output(trim_res)
                entry = _log_entry(
                    utc_now_iso(),
                    "OK" if trim_res.success else "ERR",
                    f"Trimmomatic {sample.internal_id}: {trim_res.error_message or 'OK'}" if not trim_res.success else f"Trimmomatic {sample.internal_id}: OK",
                )
                logs.append(entry)
                append_log(progress_path, entry)
                if not trim_res.success:
                    all_ok = False
                    continue
                next_sample = {
                    "internal_id": sample.internal_id,
                    "display_label": sample.display_label,
                    "mode": sample.mode.value,
                    "marker": sample.marker.value,
                    "is_control": sample.is_control,
                    "group": sample.group,
                    "r1_path": str(trim_res.files.get("r1_trimmed", sample.paths.r1)),
                    "r2_path": str(trim_res.files.get("r2_trimmed")) if trim_res.files.get("r2_trimmed") else None,
                }
                next_samples.append(next_sample)
                is_extended_frags = "extendedfrag" in sample.internal_id.lower() or "extendedfrag" in (sample.display_label or "").lower()
                if run_vsearch_branch and sample.mode.value == "PE" and not is_extended_frags:
                    if _ensure_existing_path(trim_res.files.get("r1_trimmed"), "Preprocess", sample.internal_id, "r1_trimmed") is None:
                        all_ok = False
                        continue
                    if _ensure_existing_path(trim_res.files.get("r2_trimmed"), "Preprocess", sample.internal_id, "r2_trimmed") is None:
                        all_ok = False
                        continue
                    pear_res = run_pear(sample, run_config, trim_res)
                    preprocess_outputs[sample.internal_id]["pear"] = serialize_step_output(pear_res)
                    entry = _log_entry(
                        utc_now_iso(),
                        "OK" if pear_res.success else "ERR",
                        f"PEAR {sample.internal_id}: {pear_res.error_message or 'OK'}" if not pear_res.success else f"PEAR {sample.internal_id}: OK",
                    )
                    logs.append(entry)
                    append_log(progress_path, entry)
                    if not pear_res.success:
                        all_ok = False


            entry = _log_entry(utc_now_iso(), "DBG", f"Preprocess: next_samples={[s['internal_id'] for s in next_samples]}, preprocess_outputs_keys={list(preprocess_outputs.keys())}")
            logs.append(entry)
            append_log(progress_path, entry)
            # Parse the Trimmomatic/PEAR read accounting and the merged-fragment
            # lengths once, here: the length histogram streams multi-GB FASTQs and
            # must never run inside a UI callback.
            try:
                preprocess_metrics_files = {
                    key: str(path) for key, path in build_preprocess_metrics(run_config.output_dir / "preprocess").items()
                }
            except Exception as exc:
                entry = _log_entry(utc_now_iso(), "ERR", f"Failed to build preprocess metrics: {exc}")
                logs.append(entry)
                append_log(progress_path, entry)
            preprocess_record = {
                "signature": run_config.run_id,
                "status": "success" if all_ok else "failure",
                "updated_at": utc_now_iso(),
                "outputs": preprocess_outputs,
                "next_samples": next_samples,
                "metrics_files": preprocess_metrics_files,
                "log_lines": logs.copy(),
            }
            progress_state["preprocess_record"] = preprocess_record
            stages["Preprocess"] = preprocess_record["status"]
            _write_progress_state()
            active_stage_results.append(preprocess_record["status"] == "success")
            _publish()
        else:
            preprocess_record = _skipped_record()
            progress_state["preprocess_record"] = preprocess_record
            stages["Preprocess"] = "skipped"
            _write_progress_state()
            _publish()

        if not run_preprocess:
            next_samples = [
                {
                    "internal_id": sample.internal_id,
                    "display_label": sample.display_label,
                    "mode": sample.mode.value,
                    "marker": sample.marker.value,
                    "is_control": sample.is_control,
                    "group": sample.group,
                    "r1_path": str(sample.paths.r1),
                    "r2_path": str(sample.paths.r2) if sample.paths.r2 else None,
                }
                for sample in run_config.samples
            ]

        # Analysis
        analysis_outputs = {}
        analysis_record = _skipped_record()
        # Each analysis branch contributes one feature set (OTU and/or ASV). Every
        # later stage iterates this list, so it must exist even when Analysis is off.
        feature_sets: list[dict] = []
        if run_analysis:
            stages["Analysis"] = "running"
            _publish()
            entry = _log_entry(utc_now_iso(), "RUN", "Starting Analysis stage")
            logs.append(entry)
            append_log(progress_path, entry)
            success = True
            # Which branch failed and why, kept per branch rather than collapsed into a
            # single stage-wide boolean — the stage can be half-good and still useful.
            branch_errors: dict[str, str] = {}
            # Samples that fell out of the pool. Surfaced loudly, but not fatal.
            dropped_samples: list[str] = []
            if run_vsearch_branch:
                filtered_fastas: dict[str, Path] = {}
                for sample_record in next_samples:
                    sample_id = sample_record["internal_id"]
                    sample_mode = str(sample_record.get("mode") or "SE").upper()
                    if run_preprocess and sample_mode == "PE":
                        pear_record = (preprocess_outputs.get(sample_id) or {}).get("pear")
                        if not pear_record:
                            assembled_value = sample_record.get("r1_path")
                        else:
                            assembled_value = (pear_record.get("files") or {}).get("assembled")
                    else:
                        assembled_value = sample_record.get("r1_path")

                    # A sample that cannot be filtered drops out of the pool. That is
                    # worth shouting about — its reads are simply absent from the OTU
                    # table — but it must not discard the samples that did work.
                    if not assembled_value:
                        dropped_samples.append(f"{sample_id} (no input)")
                        entry = _log_entry(utc_now_iso(), "ERR", f"No input FASTA/FASTQ for {sample_id}")
                        logs.append(entry)
                        append_log(progress_path, entry)
                        continue

                    assembled_path = Path(assembled_value)
                    if not assembled_path.exists():
                        dropped_samples.append(f"{sample_id} (input does not exist)")
                        entry = _log_entry(utc_now_iso(), "ERR", f"Input missing: {assembled_path}")
                        logs.append(entry)
                        append_log(progress_path, entry)
                        continue

                    filter_res = run_vsearch_filter_sample(
                        assembled_path,
                        run_config,
                        run_config.output_dir / "analysis" / sample_id,
                        sample_id,
                    )
                    analysis_outputs[sample_id] = {"vsearch_filter": serialize_step_output(filter_res)}
                    if filter_res.success:
                        filtered_fastas[sample_id] = filter_res.files["filtered_fasta"]
                        entry = _log_entry(utc_now_iso(), "OK", f"VSEARCH filter {sample_id}: OK")
                    else:
                        dropped_samples.append(f"{sample_id} (filtro)")
                        entry = _log_entry(utc_now_iso(), "ERR", f"VSEARCH filter {sample_id}: {filter_res.error_message or 'unknown error'}")
                    logs.append(entry)
                    append_log(progress_path, entry)

                if filtered_fastas:
                    pool_res = run_vsearch_pool_and_cluster(filtered_fastas, run_config, run_config.output_dir / "analysis" / "pooled")
                    analysis_outputs["vsearch_pool"] = serialize_step_output(pool_res)
                    if pool_res.success:
                        combined_fasta = pool_res.files["otus_fasta"]
                        analysis_outputs["combined_fasta_vsearch"] = str(combined_fasta)
                        entry = _log_entry(utc_now_iso(), "OK", "VSEARCH pooled clustering: OK")
                        logs.append(entry)
                        append_log(progress_path, entry)
                        otu_table_tsv = pool_res.files.get("otu_table_tsv")
                        if otu_table_tsv and Path(otu_table_tsv).exists():
                            analysis_outputs["otu_table_tsv"] = str(otu_table_tsv)
                            analysis_outputs["combined_otu_uc"] = str(pool_res.files.get("map_uc", ""))
                            try:
                                matrix = pd.read_csv(otu_table_tsv, sep="\t")
                                matrix_path = run_config.output_dir / "analysis" / "otu_matrix.csv"
                                matrix.to_csv(matrix_path, index=False)
                                analysis_outputs["otu_matrix"] = str(matrix_path)
                            except Exception as exc:
                                entry = _log_entry(utc_now_iso(), "ERR", f"Failed to convert otu_table to CSV: {exc}")
                                logs.append(entry)
                                append_log(progress_path, entry)
                    else:
                        branch_errors["otu"] = pool_res.error_message or "unknown error"
                        analysis_outputs["combined_fasta_vsearch"] = ""
                        entry = _log_entry(utc_now_iso(), "ERR", f"VSEARCH pooled clustering: {pool_res.error_message or 'unknown error'}")
                        logs.append(entry)
                        append_log(progress_path, entry)
                else:
                    branch_errors["otu"] = "no sample produced filtered reads"
                    analysis_outputs["combined_fasta_vsearch"] = ""
                    entry = _log_entry(utc_now_iso(), "ERR", "Analysis: no sample produced filtered reads for VSEARCH.")
                    logs.append(entry)
                    append_log(progress_path, entry)

            if run_dada2_branch:
                dada2_next_samples = []
                for sample_record in next_samples:
                    trim_record = (preprocess_outputs.get(sample_record["internal_id"]) or {}).get("trimmomatic") or {}
                    trim_files = trim_record.get("files") or {}
                    r1_trim = trim_files.get("r1_trimmed")
                    r2_trim = trim_files.get("r2_trimmed")
                    if not r1_trim or not Path(r1_trim).exists():
                        r1_trim = sample_record.get("r1_path")
                        r2_trim = None
                    dada2_record = dict(sample_record)
                    dada2_record["r1_path"] = r1_trim
                    dada2_record["r2_path"] = r2_trim if r2_trim and Path(str(r2_trim)).exists() else None
                    dada2_next_samples.append(dada2_record)

                analysis_samples = [build_sample_from_record(sample_record) for sample_record in dada2_next_samples]
                dada2_result = run_dada2(analysis_samples, run_config, run_config.output_dir / "analysis" / "dada2")
                analysis_outputs["dada2"] = serialize_step_output(dada2_result)
                if dada2_result.success:
                    analysis_outputs["combined_fasta_dada2"] = str(dada2_result.files.get("asv_sequences_fasta"))
                    entry = _log_entry(utc_now_iso(), "OK", "DADA2 completed")
                    logs.append(entry)
                    append_log(progress_path, entry)
                else:
                    # One branch dying must not take the other down with it. VSEARCH can
                    # have produced a perfectly good OTU set, and throwing that away
                    # because DADA2 ran out of memory loses a run that half-worked.
                    branch_errors["asv"] = dada2_result.error_message or "unknown error"
                    entry = _log_entry(utc_now_iso(), "ERR", f"DADA2 failed: {dada2_result.error_message or 'unknown error'}")
                    logs.append(entry)
                    append_log(progress_path, entry)

            if run_vsearch_branch:
                vsearch_fasta = analysis_outputs.get("combined_fasta_vsearch") or None
                if vsearch_fasta:
                    vsearch_path = Path(vsearch_fasta)
                    vsearch_fasta = str(vsearch_path) if vsearch_path.exists() else None
                analysis_outputs["combined_fasta_vsearch"] = vsearch_fasta or ""
            if not run_vsearch_branch and not run_dada2_branch:
                success = False
                entry = _log_entry(utc_now_iso(), "ERR", f"Invalid analysis method: {analysis_method}")
                logs.append(entry)
                append_log(progress_path, entry)

            # Each analysis branch yields its own feature set. They are carried
            # forward SEPARATELY — the same taxonomy, tree and metrics are built once
            # per set — because an OTU (97% cluster) and an ASV (exact variant) are
            # different units and must never share a table.
            feature_sets = feature_sets_from_outputs(analysis_outputs)
            analysis_outputs["feature_sets"] = [
                {"key": fs["key"], "label": fs["label"], "fasta": str(fs["fasta"]), "matrix": str(fs["matrix"]) if fs["matrix"] else ""}
                for fs in feature_sets
            ]
            # Legacy single-track keys still resolve, to the primary (first) set.
            analysis_outputs["combined_fasta"] = str(feature_sets[0]["fasta"]) if feature_sets else ""
            for fs in feature_sets:
                entry = _log_entry(
                    utc_now_iso(),
                    "OK",
                    f"Feature set {fs['label']}: {_count_features(fs['fasta'])} features -> {fs['key']}/",
                )
                logs.append(entry)
                append_log(progress_path, entry)

            # The stage succeeds when it produced SOMETHING to classify. A branch that
            # died is a degradation, not a write-off: losing a good VSEARCH OTU set
            # because DADA2 hit the OOM killer means the run is discarded even though
            # half of it worked, which is the worst possible outcome for the user.
            analysis_outputs["branch_errors"] = dict(branch_errors)
            analysis_outputs["dropped_samples"] = list(dropped_samples)
            if dropped_samples:
                entry = _log_entry(
                    utc_now_iso(),
                    "WARN",
                    f"Samples excluded from the feature table: {', '.join(dropped_samples)}. "
                    f"The reads from these samples are NOT in the results.",
                )
                logs.append(entry)
                append_log(progress_path, entry)
            if feature_sets and branch_errors:
                degraded = ", ".join(f"{key.upper()}: {message}" for key, message in branch_errors.items())
                entry = _log_entry(
                    utc_now_iso(),
                    "WARN",
                    f"Partial analysis — continuing with {len(feature_sets)} set(s). Failed → {degraded}",
                )
                logs.append(entry)
                append_log(progress_path, entry)
            success = success and bool(feature_sets)

            # Rarefaction is a Monte-Carlo subsampling over a grid of depths; it has
            # to be computed here, once, and never inside a render callback.
            analysis_metrics_files: dict[str, str] = {}
            try:
                analysis_metrics_files = build_metrics_from_outputs(run_config.output_dir / "analysis", analysis_outputs)
            except Exception as exc:
                entry = _log_entry(utc_now_iso(), "ERR", f"Failed to build analysis metrics: {exc}")
                logs.append(entry)
                append_log(progress_path, entry)

            # Rebuild the read-tracking waterfall now that the Analysis stage exists, so
            # it can close the account: the quality filter and the OTU mapping both drop
            # reads, and until they were included the table stopped at "assembled" and
            # ~22k reads per run simply went unmentioned. Same builder as the Preprocess
            # stage, given the analysis dir this time.
            try:
                refreshed = {
                    key: str(path)
                    for key, path in build_preprocess_metrics(
                        run_config.output_dir / "preprocess", run_config.output_dir / "analysis"
                    ).items()
                }
                if refreshed:
                    preprocess_metrics_files.update(refreshed)
                    preprocess_record["metrics_files"] = dict(preprocess_metrics_files)
                    progress_state["preprocess_record"] = preprocess_record
            except Exception as exc:
                entry = _log_entry(utc_now_iso(), "ERR", f"Failed to complete read tracking with the analysis: {exc}")
                logs.append(entry)
                append_log(progress_path, entry)
            analysis_record = {
                "signature": run_config.run_id,
                "status": "success" if success else "failure",
                "updated_at": utc_now_iso(),
                "outputs": analysis_outputs,
                "metrics_files": analysis_metrics_files,
                "branch_errors": dict(branch_errors),
                "log_lines": logs.copy(),
            }
            progress_state["analysis_record"] = analysis_record
            stages["Analysis"] = analysis_record["status"]
            _write_progress_state()
            active_stage_results.append(analysis_record["status"] == "success")
            _publish()
        else:
            stages["Analysis"] = "skipped"
            progress_state["analysis_record"] = analysis_record
            _write_progress_state()
            _publish()

        # Databases — run once per feature set, into its own directory.
        databases_outputs: dict = {}
        db_record = _skipped_record()
        if run_databases:
            stages["Databases"] = "running"
            _publish()
            if not feature_sets:
                entry = _log_entry(utc_now_iso(), "ERR", "Databases: the analysis produced no feature set (OTU/ASV)")
                logs.append(entry)
                append_log(progress_path, entry)
                db_record = _failed_record(run_config.run_id, databases_outputs, logs)
                stages["Databases"] = "failure"
                active_stage_results.append(False)
                progress_state["db_record"] = db_record
                _write_progress_state()
                _publish()
            elif analysis_record.get("status") not in ("success", "skipped"):
                entry = _log_entry(utc_now_iso(), "ERR", "Databases: the previous analysis failed")
                logs.append(entry)
                append_log(progress_path, entry)
                db_record = _failed_record(run_config.run_id, databases_outputs, logs)
                stages["Databases"] = "failure"
                active_stage_results.append(False)
                progress_state["db_record"] = db_record
                _write_progress_state()
                _publish()
            else:
                entry = _log_entry(utc_now_iso(), "RUN", "Starting Databases stage")
                logs.append(entry)
                append_log(progress_path, entry)
                db_sources = {
                    str(source).lower()
                    for source in (
                        (run_config_data or {}).get("db_sources")
                        or run_config.vsearch_params.get("db_sources")
                        or []
                    )
                }
                db_paths = get_db_paths()
                sintax_cutoff = float(run_config.vsearch_params.get("sintax_cutoff", 0.8))
                custom_db_raw = (run_config.vsearch_params or {}).get("custom_db_path")
                custom_db_path = Path(custom_db_raw) if custom_db_raw and Path(custom_db_raw).exists() else None
                custom_db_kind = (run_config.vsearch_params or {}).get("custom_db_kind")

                # The marker and the primers this run amplified with are what SILVA/PR2
                # get cut down to, so the reference describes the same stretch of gene as
                # the queries — and fits in RAM on the way.
                run_marker = _pipeline_run_marker(run_config)
                run_primers = resolve_primers(run_config)

                # Report up front what will not fit in RAM, instead of letting the OS
                # kill a source halfway through and calling that a "failure".
                for note in preflight_notes(
                    db_sources,
                    db_paths,
                    len(run_config.samples),
                    marker_filter=run_marker if marker_filter_applies(run_marker, run_primers) else "",
                ):
                    entry = _log_entry(utc_now_iso(), "WARN", note)
                    logs.append(entry)
                    append_log(progress_path, entry)

                silva_config = SILVAConfig(db_path=db_paths.silva) if "silva" in db_sources and db_paths.silva and db_paths.silva.exists() else None
                pr2_config = PR2Config(db_path=db_paths.pr2, tax_path=db_paths.pr2_tax) if "pr2" in db_sources and db_paths.pr2 and db_paths.pr2.exists() else None
                # NCBI's public remote nt is the default source, by decision: it is what a
                # reviewer can reproduce without our filesystem. A local db is an opt-in
                # emergency exit (PARAFISH_NCBI_LOCAL_DB) for when NCBI is down — it used
                # to be picked automatically whenever the files happened to be present,
                # which silently changed the source out from under the run.
                #
                # The task differs by mode on purpose. The public quota only completes
                # megablast (blastn/dc-megablast are refused with SIGXCPU); a local db is
                # our own CPU and can afford the sensitive search that finds a parasite at
                # ~90% identity, which megablast's 28 bp seed can miss.
                use_local_blast = NCBI_USE_LOCAL_DB and bool(db_paths.blast and db_paths.blast.exists())
                ncbi_config = (
                    NCBIConfig(
                        local_db_path=db_paths.blast if use_local_blast else None,
                        use_remote=not use_local_blast,
                        task=str(run_config.vsearch_params.get("ncbi_task") or ("blastn" if use_local_blast else "megablast")),
                        word_size=int(run_config.vsearch_params["ncbi_word_size"]) if run_config.vsearch_params.get("ncbi_word_size") else (11 if use_local_blast else None),
                        evalue=float(run_config.vsearch_params.get("ncbi_evalue", 1e-10)),
                        max_target_seqs=int(run_config.vsearch_params.get("ncbi_maxhits", 20)),
                    )
                    if "ncbi" in db_sources
                    else None
                )
                if "ncbi" in db_sources:
                    entry = _log_entry(
                        utc_now_iso(),
                        "RUN",
                        f"NCBI: {'local database' if use_local_blast else 'public remote (nt)'} · task={ncbi_config.task}",
                    )
                    logs.append(entry)
                    append_log(progress_path, entry)
                # BOLDigger3 uses the remote BOLD API; do not gate on a local dir.
                bold_config = (
                    BOLDConfig(
                        db_path=db_paths.bold_dir or Path("databases/bold"),
                        db_source=int(run_config.bold_params.get("db", 1)),
                        mode=int(run_config.bold_params.get("mode", 2)),
                    )
                    if "bold" in db_sources
                    else None
                )

                taxonomy_by_set: dict = {}
                db_metrics_by_set: dict = {}
                any_success = False
                for feature_set in feature_sets:
                    key = feature_set["key"]
                    fasta_path = feature_set["fasta"]
                    abundance_matrix_path = feature_set["matrix"]
                    entry = _log_entry(utc_now_iso(), "RUN", f"Taxonomy for {feature_set['label']}")
                    logs.append(entry)
                    append_log(progress_path, entry)

                    tax_result = run_taxonomy_pipeline(
                        fasta_path,
                        run_config.output_dir / "databases" / key,
                        silva_config=silva_config,
                        pr2_config=pr2_config,
                        ncbi_config=ncbi_config,
                        bold_config=bold_config,
                        custom_db_path=custom_db_path,
                        custom_db_kind=custom_db_kind,
                        abundance_matrix_path=abundance_matrix_path,
                        consensus_method=str(run_config.vsearch_params.get("db_consensus") or "lca"),
                        sintax_cutoff=sintax_cutoff,
                        threads=int(run_config.vsearch_params.get("threads", 8)),
                        feature_set=key,
                        marker=run_marker,
                        primers=run_primers,
                    )
                    taxonomy_by_set[key] = serialize_step_output(tax_result)
                    any_success = any_success or tax_result.success
                    entry = _log_entry(
                        utc_now_iso(),
                        "OK" if tax_result.success else "ERR",
                        f"Taxonomy {key}: {tax_result.error_message or 'OK'}" if not tax_result.success else f"Taxonomy {key} completed",
                    )
                    logs.append(entry)
                    append_log(progress_path, entry)

                    # Resolution, decontamination and prevalence are computed here so
                    # the UI never re-runs PERMANOVA or Fisher tests on a page switch.
                    try:
                        db_metrics_by_set[key] = {
                            name: str(path)
                            for name, path in build_database_metrics(
                                run_config.output_dir / "databases" / key,
                                tax_result.files.get("taxonomy_consensus_csv"),
                                abundance_matrix_path,
                                control_samples=[sample.internal_id for sample in run_config.samples if sample.is_control],
                                groups={sample.internal_id: sample.group for sample in run_config.samples if sample.group},
                            ).items()
                        }
                    except Exception as exc:
                        entry = _log_entry(utc_now_iso(), "ERR", f"Failed to build taxonomy metrics ({key}): {exc}")
                        logs.append(entry)
                        append_log(progress_path, entry)

                databases_outputs["taxonomy_by_set"] = taxonomy_by_set

                # OTU vs ASV, side by side. Two separate summaries plus a bridge that
                # joins on the TAXON, never on the feature id — that is the only
                # comparison between the two that is biologically meaningful.
                if len(feature_sets) > 1:
                    try:
                        consensus_by_set = {}
                        for key, output in taxonomy_by_set.items():
                            consensus = (output.get("files") or {}).get("taxonomy_consensus_csv")
                            if consensus:
                                consensus_by_set[key] = Path(consensus)
                        comparison_files = build_feature_set_comparison(
                            run_config.output_dir / "analysis" / "comparison",
                            feature_sets,
                            consensus_by_set=consensus_by_set,
                        )
                        databases_outputs["feature_set_comparison"] = {k: str(v) for k, v in comparison_files.items()}
                        entry = _log_entry(utc_now_iso(), "OK", "OTU vs ASV comparison written (separate tables + taxon bridge)")
                        logs.append(entry)
                        append_log(progress_path, entry)
                    except Exception as exc:
                        entry = _log_entry(utc_now_iso(), "ERR", f"Failed to compare OTU/ASV: {exc}")
                        logs.append(entry)
                        append_log(progress_path, entry)

                # Legacy single-track key resolves to the primary feature set.
                primary_key = feature_sets[0]["key"]
                databases_outputs["taxonomy"] = taxonomy_by_set.get(primary_key, {})
                db_metrics_files = dict(db_metrics_by_set.get(primary_key, {}))
                for key, metrics in db_metrics_by_set.items():
                    for name, path in metrics.items():
                        db_metrics_files[f"{key}.{name}"] = path

                db_record = {
                    "signature": run_config.run_id,
                    "status": "success" if any_success else "failure",
                    "updated_at": utc_now_iso(),
                    "outputs": databases_outputs,
                    "metrics_files": db_metrics_files,
                    "log_lines": logs.copy(),
                }
                progress_state["db_record"] = db_record
                stages["Databases"] = db_record["status"]
                _write_progress_state()
                active_stage_results.append(db_record["status"] == "success")
                _publish()
        else:
            stages["Databases"] = "skipped"
            progress_state["db_record"] = db_record
            _write_progress_state()
            _publish()

        # Phylo
        phylo_outputs = {}
        phylo_record = _skipped_record()
        if run_phylo_stage:
            stages["Phylo"] = "running"
            _publish()
            # One tree per feature set: an OTU tree and an ASV tree are different
            # trees over different tips and cannot be superimposed.
            taxonomy_by_set = (db_record.get("outputs") or {}).get("taxonomy_by_set") or {}
            phylo_by_set: dict = {}
            phylo_metrics_by_set: dict = {}
            any_phylo_success = False
            phylo_attempted = False

            for feature_set in feature_sets:
                key = feature_set["key"]
                # The annotated FASTA carries the taxonomy in its headers, so tips are
                # readable; fall back to the bare representatives when Databases is off.
                annotated = ((taxonomy_by_set.get(key) or {}).get("files") or {}).get("annotated_fasta")
                phylo_path = Path(annotated) if annotated and Path(annotated).exists() else feature_set["fasta"]

                try:
                    seq_count = sum(1 for _ in SeqIO.parse(str(phylo_path), "fasta"))
                except Exception:
                    seq_count = 0
                guard_error = phylo_input_error(seq_count)
                if guard_error:
                    entry = _log_entry(utc_now_iso(), "ERR", f"Phylo {key}: {guard_error}")
                    logs.append(entry)
                    append_log(progress_path, entry)
                    continue

                phylo_attempted = True
                entry = _log_entry(utc_now_iso(), "RUN", f"Tree for {feature_set['label']}")
                logs.append(entry)
                append_log(progress_path, entry)
                phylo_result = run_phylo_step(phylo_path, run_config, run_config.output_dir / "phylo" / key)
                phylo_by_set[key] = {
                    "files": {k: str(v) for k, v in getattr(phylo_result, "files", {}).items()},
                    "metrics": getattr(phylo_result, "metrics", {}),
                }
                any_phylo_success = any_phylo_success or phylo_result.success
                entry = _log_entry(
                    utc_now_iso(),
                    "OK" if phylo_result.success else "ERR",
                    f"Phylo {key}: {phylo_result.error_message or 'OK'}" if not phylo_result.success else f"Phylo {key} completed",
                )
                logs.append(entry)
                append_log(progress_path, entry)

                # UniFrac is O(tips x samples) per pair; compute it here, once, and
                # against this feature set's own abundance matrix.
                try:
                    phylo_metrics_by_set[key] = build_metrics_from_files(
                        run_config.output_dir / "phylo" / key,
                        phylo_by_set[key]["files"],
                        str(feature_set["matrix"]) if feature_set["matrix"] else None,
                    )
                except Exception as exc:
                    entry = _log_entry(utc_now_iso(), "ERR", f"Failed to build phylogeny metrics ({key}): {exc}")
                    logs.append(entry)
                    append_log(progress_path, entry)

            if not phylo_attempted:
                entry = _log_entry(utc_now_iso(), "ERR", "Phylo: no feature set has enough sequences")
                logs.append(entry)
                append_log(progress_path, entry)
                phylo_record = _failed_record(run_config.run_id, phylo_outputs, logs)
                stages["Phylo"] = "failure"
                active_stage_results.append(False)
                progress_state["phylo_record"] = phylo_record
                _write_progress_state()
                _publish()
            else:
                phylo_outputs["phylo_by_set"] = phylo_by_set
                primary_key = feature_sets[0]["key"]
                phylo_outputs["phylo"] = phylo_by_set.get(primary_key, {})
                phylo_metrics_files = dict(phylo_metrics_by_set.get(primary_key, {}))
                for key, metrics in phylo_metrics_by_set.items():
                    for name, path in metrics.items():
                        phylo_metrics_files[f"{key}.{name}"] = path
                phylo_record = {
                    "signature": run_config.run_id,
                    "status": "success" if any_phylo_success else "failure",
                    "updated_at": utc_now_iso(),
                    "outputs": phylo_outputs,
                    "metrics_files": phylo_metrics_files,
                    "log_lines": logs.copy(),
                }
                progress_state["phylo_record"] = phylo_record
                stages["Phylo"] = phylo_record["status"]
                _write_progress_state()
                active_stage_results.append(phylo_record["status"] == "success")
                _publish()
        else:
            stages["Phylo"] = "skipped"
            progress_state["phylo_record"] = phylo_record
            _write_progress_state()
            _publish()

        overall = "success" if all(active_stage_results) else "failure"
        final_entry = _log_entry(utc_now_iso(), "OK" if overall == "success" else "ERR", f"Pipeline finished: {overall}")
        logs.append(final_entry)
        append_log(progress_path, final_entry)
        progress_state.update(
            {
                "done": True,
                "overall": overall,
                "stages": dict(stages),
                "qc_record": qc_record,
                "preprocess_record": preprocess_record,
                "analysis_record": analysis_record,
                "db_record": db_record,
                "phylo_record": phylo_record,
            }
        )
        write_progress(progress_path, progress_state)
    except Exception as exc:
        entry = _log_entry(utc_now_iso(), "ERR", f"Pipeline crashed: {exc}")
        logs.append(entry)
        append_log(progress_path, entry)
        stages["Pipeline"] = "failure"
        progress_state.update(
            {
                "done": True,
                "overall": "failure",
                "stages": dict(stages),
                "qc_record": qc_record,
                "preprocess_record": preprocess_record,
                "analysis_record": analysis_record,
                "db_record": db_record,
                "phylo_record": phylo_record,
            }
        )
        write_progress(progress_path, progress_state)
    finally:
        try:
            with _progress_lock:
                _pipeline_threads.pop(run_config.run_id, None)
            pdir = str(run_config.output_dir.parent.parent)
            _apply_cleanup_post_run(pdir, str(run_config.output_dir))
        except Exception:
            pass


def _any_pipeline_running() -> bool:
    """True iff a pipeline thread is genuinely alive in THIS process.

    The single source of truth for "is a run happening" — the live thread, never
    progress.json. An OOM-killed run leaves ``stages: running`` (and no ``done``) on
    disk forever; a disk-based check would then refuse every future run and disable
    the Start button permanently. The process alone knows what is really executing.
    """
    with _progress_lock:
        return any(thread.is_alive() for thread in _pipeline_threads.values())


def _launch_pipeline_run(run_config, run_config_data) -> bool:
    """Starts the background pipeline thread for run_config, unless a run is already
    live in this process.

    The sole launcher, wired only to the Progress page's explicit "Start Pipeline"
    button. The guard is atomic (check-and-start under the lock) so a fast double
    click or a second browser tab cannot slip a second thread past it, and it refuses
    ANY concurrent run — not just the same run_id — because on this box two pipelines
    at once is a guaranteed OOM. It keys on the live thread, never on progress.json,
    so a stale "running" left by an OS-killed run does not block a fresh launch.

    Returns True if a new thread was started, False if one was already running.
    """
    progress_path = run_config.output_dir / "progress.json"
    with _progress_lock:
        if any(thread.is_alive() for thread in _pipeline_threads.values()):
            return False
        thread = threading.Thread(
            target=_run_pipeline_thread,
            args=(run_config, run_config_data, progress_path),
            daemon=True,
        )
        _pipeline_threads[run_config.run_id] = thread
        thread.start()
    return True


def register_callbacks(app: Any) -> None:
    @app.callback(
        Output("progress-interval", "disabled", allow_duplicate=True),
        Output("progress-overall-status", "children", allow_duplicate=True),
        Output("btn-start-progress", "disabled", allow_duplicate=True),
        Input("btn-start-progress", "n_clicks"),
        State("store-project", "data"),
        State("store-samples", "data"),
        State("store-run-config", "data"),
        prevent_initial_call=True,
    )
    def _start_pipeline(n_clicks, project_data, samples_data, run_config_data):
        if not n_clicks:
            raise PreventUpdate
        run_config = run_config_from_store(project_data, samples_data, run_config_data)
        run_config.samples = _selected_samples(run_config.samples)
        started = _launch_pipeline_run(run_config, run_config_data)
        # Disable the button the instant a run is live, so a second click or tab cannot
        # even attempt a launch; the backend guard is the real lock, this is the UI half.
        return False, "Running..." if started else "Already running...", _any_pipeline_running()

    @app.callback(
        Output("progress-stage-rows", "children"),
        Output("progress-logs", "children"),
        Output("progress-overall-status", "children", allow_duplicate=True),
        Output("progress-interval", "disabled", allow_duplicate=True),
        Output("btn-start-progress", "disabled", allow_duplicate=True),
        Output("store-qc-outputs", "data", allow_duplicate=True),
        Output("store-preprocess-outputs", "data", allow_duplicate=True),
        Output("store-analysis-outputs", "data", allow_duplicate=True),
        Output("store-databases-outputs", "data", allow_duplicate=True),
        Output("store-phylo-outputs", "data", allow_duplicate=True),
        Input("progress-interval", "n_intervals"),
        State("store-run-config", "data"),
        State("store-project", "data"),
        State("store-samples", "data"),
        State("store-qc-outputs", "data"),
        State("store-preprocess-outputs", "data"),
        State("store-analysis-outputs", "data"),
        State("store-databases-outputs", "data"),
        State("store-phylo-outputs", "data"),
        prevent_initial_call=True,
    )
    def _poll_progress(
        n_intervals,
        run_config_data,
        project_data,
        samples_data,
        qc_current,
        preprocess_current,
        analysis_current,
        databases_current,
        phylo_current,
    ):
        run_id = (run_config_data or {}).get("run_id")
        if not run_id:
            raise PreventUpdate

        run_config = run_config_from_store(project_data, samples_data, run_config_data)
        run_config.samples = _selected_samples(run_config.samples)
        progress_path = run_config.output_dir / "progress.json"
        state = read_progress(progress_path)
        if not state:
            raise PreventUpdate
        if state.get("run_id") != run_id:
            raise PreventUpdate

        all_log_lines = state.get("log_lines") or []
        seen = set()
        log_lines = []
        for entry in all_log_lines:
            key = (entry.get("ts"), entry.get("level"), entry.get("msg"))
            if key not in seen:
                seen.add(key)
                log_lines.append(entry)
        stages = state.get("stages") or {}
        done = state.get("done", False)
        overall = state.get("overall", "running")

        stage_rows = [html.Div(_stage_entry(name, status)) for name, status in stages.items()]
        log_block = html.Pre("\n".join(f"[{line.get('ts')}] {line.get('level')}: {line.get('msg')}" for line in log_lines))
        overall_text = f"Overall: {overall}" if done else "Running..."

        # The stage pages (QC, Preprocess, Analysis, Databases, Phylo) fill in only once
        # the WHOLE run has finished — not one at a time as each stage completes. During
        # the run the Progress page shows the live log and the stage-status rows; the
        # heavy per-stage tables and plots appear together at the end.
        #
        # This is crash-safe: the pipeline thread's `except` handler also sets
        # ``done: True`` and writes every stage record, so a run that fails still
        # publishes whatever succeeded — you are not left staring at empty pages after a
        # crash. It is only a hard OS kill of the whole process (OOM) that would leave
        # nothing, and for that the Results page reads the run straight from disk.
        # After `done`, the interval keeps ticking for a moment before it disables. Only
        # publish a record the first time it differs from what the store already holds,
        # so those trailing ticks do not re-render the heavy stage tables and plots.
        def _first_change(record, current):
            # Stamped with the run that produced it, so a stage page can tell this
            # record from the one the PREVIOUS run left in the same session store.
            stamped = stamp_record(record, run_id)
            return stamped if stamped != current else no_update

        # The Start button follows the LIVE THREAD, not this disk snapshot: it stays
        # disabled while a run is genuinely executing and re-enables the moment the
        # thread ends — including an OOM kill, which never writes `done` and would
        # otherwise leave the button dead forever if it keyed on progress.json.
        start_disabled = _any_pipeline_running()

        if not done:
            return stage_rows, log_block, overall_text, False, start_disabled, no_update, no_update, no_update, no_update, no_update

        return (
            stage_rows,
            log_block,
            overall_text,
            True,
            start_disabled,
            _first_change(state.get("qc_record"), qc_current),
            _first_change(state.get("preprocess_record"), preprocess_current),
            _first_change(state.get("analysis_record"), analysis_current),
            _first_change(state.get("db_record"), databases_current),
            _first_change(state.get("phylo_record"), phylo_current),
        )