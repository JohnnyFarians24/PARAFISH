"""Orchestration callbacks for running the full pipeline from the Setup page."""

from __future__ import annotations

import json
import threading
import subprocess
from pathlib import Path
from typing import Any

from Bio import SeqIO
from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from app import _pipeline_threads, _progress_lock
from analysis.vsearch_runner import _parse_vsearch_uc_to_matrix
from analysis.dada2_runner import run_dada2
from analysis.vsearch_runner import run_vsearch
from core.config import get_db_paths, VSEARCH_PATH
from core.progress_writer import append_log, read_progress, write_progress
from core.storage import _apply_cleanup_post_run
from core.utils import prefix_fasta_ids
from databases.bold_db import BOLDConfig
from databases.ncbi_db import NCBIConfig
from databases.pr2_db import PR2Config
from databases.silva_db import SILVAConfig
from databases.taxonomy_runner import run_taxonomy_pipeline
from phylo.phylo_runner import run_phylo as run_phylo_step
from preprocess.pear_runner import run_pear
from preprocess.trimmomatic_runner import run_trimmomatic
from qc.fastqc_runner import run_fastqc
from ui.pipeline_state import utc_now_iso
from ui.runtime_state import build_sample_from_record, resolve_default_sintax_database, run_config_from_store, serialize_step_output

def _stage_entry(stage_name: str, status: str) -> str:
    return f"{stage_name}: {status}"


def _log_entry(ts: str, level: str, msg: str) -> dict:
    return {"ts": ts, "level": level, "msg": msg}


def _selected_samples(samples):
	return [sample for sample in samples or [] if getattr(sample, "selected", True)]


def _run_pipeline_thread(run_config, run_config_data, progress_path):
    """Runs in a background thread and writes progress to disk."""
    del run_config_data

    stages: dict[str, str] = {}
    logs: list[dict] = []
    active_stage_results: list[bool] = []
    progress_state = {
        "run_id": run_config.run_id,
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

    def _publish() -> None:
        return None

    def _write_progress_state() -> None:
        progress_state.update({
            "stages": dict(stages),
            "log_lines": logs.copy(),
        })
        write_progress(progress_path, progress_state)

    def _skipped_record() -> dict:
        return {"status": "skipped", "outputs": {}, "updated_at": utc_now_iso()}

    def _ensure_existing_path(path_value, stage_name: str, sample_id: str, label: str):
        if path_value is None or path_value == "":
            entry = _log_entry(utc_now_iso(), "ERR", f"{stage_name}: {label} ausente para {sample_id}")
            logs.append(entry)
            append_log(progress_path, entry)
            return None
        candidate = Path(path_value)
        if not candidate.exists():
            entry = _log_entry(utc_now_iso(), "ERR", f"{stage_name}: {label} inexistente para {sample_id}: {candidate}")
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
            qc_record = {
                "signature": run_config.run_id,
                "status": "success" if all_ok else "failure",
                "updated_at": utc_now_iso(),
                "outputs": qc_outputs,
                "log_lines": logs.copy(),
            }
            progress_state["qc_record"] = qc_record
            _write_progress_state()
            stages["QC"] = qc_record["status"]
            active_stage_results.append(qc_record["status"] == "success")
            _publish()
        else:
            qc_record = _skipped_record()
            progress_state["qc_record"] = qc_record
            _write_progress_state()
            stages["QC"] = "skipped"
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
                trim_res = run_trimmomatic(sample, run_config)
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
            preprocess_record = {
                "signature": run_config.run_id,
                "status": "success" if all_ok else "failure",
                "updated_at": utc_now_iso(),
                "outputs": preprocess_outputs,
                "next_samples": next_samples,
                "log_lines": logs.copy(),
            }
            progress_state["preprocess_record"] = preprocess_record
            _write_progress_state()
            stages["Preprocess"] = preprocess_record["status"]
            active_stage_results.append(preprocess_record["status"] == "success")
            _publish()
        else:
            preprocess_record = _skipped_record()
            progress_state["preprocess_record"] = preprocess_record
            _write_progress_state()
            stages["Preprocess"] = "skipped"
            _publish()

        if not run_preprocess:
            next_samples = [
                {
                    "internal_id": sample.internal_id,
                    "display_label": sample.display_label,
                    "mode": sample.mode.value,
                    "marker": sample.marker.value,
                    "r1_path": str(sample.paths.r1),
                    "r2_path": str(sample.paths.r2) if sample.paths.r2 else None,
                }
                for sample in run_config.samples
            ]

        # Analysis
        analysis_outputs = {}
        analysis_record = _skipped_record()
        if run_analysis:
            stages["Analysis"] = "running"
            _publish()
            entry = _log_entry(utc_now_iso(), "RUN", "Starting Analysis stage")
            logs.append(entry)
            append_log(progress_path, entry)
            success = True
            if run_vsearch_branch:
                combined_fasta = run_config.output_dir / "analysis" / "combined_otus.fasta"
                combined_fasta.parent.mkdir(parents=True, exist_ok=True)
                with combined_fasta.open("w", encoding="utf-8") as combined_handle:
                    for sample_record in next_samples:
                        sample_id = sample_record["internal_id"]
                        sample_mode = str(sample_record.get("mode") or "SE").upper()
                        sample_marker = str(sample_record.get("marker") or "")
                        if run_preprocess and sample_mode == "PE":
                            pear_record = (preprocess_outputs.get(sample_id) or {}).get("pear")
                            if not pear_record:
                                assembled_value = sample_record.get("r1_path")
                            else:
                                assembled_value = (pear_record.get("files") or {}).get("assembled")
                        else:
                            assembled_value = sample_record.get("r1_path")

                        if not assembled_value:
                            success = False
                            entry = _log_entry(utc_now_iso(), "ERR", f"No input FASTA/FASTQ for {sample_id}")
                            logs.append(entry)
                            append_log(progress_path, entry)
                            continue

                        assembled_path = Path(assembled_value)
                        if not assembled_path.exists():
                            success = False
                            entry = _log_entry(utc_now_iso(), "ERR", f"Input missing: {assembled_path}")
                            logs.append(entry)
                            append_log(progress_path, entry)
                            continue
        

                        db_sintax_path = resolve_default_sintax_database(sample_marker)
                        entry = _log_entry(utc_now_iso(), "DBG", f"VSEARCH input para {sample_id}: {assembled_path} (existe={assembled_path.exists()}, bytes={assembled_path.stat().st_size if assembled_path.exists() else 'N/A'})")
                        logs.append(entry)
                        append_log(progress_path, entry)
                        vsearch_res = run_vsearch(
                            assembled_path,
                            run_config,
                            run_config.output_dir / "analysis" / sample_id,
                            db_sintax=Path(db_sintax_path) if db_sintax_path else None,
                        )
                        analysis_outputs[sample_id] = {"vsearch": serialize_step_output(vsearch_res)}
                        if vsearch_res.success:
                            otu_path = vsearch_res.files.get("otus_fasta")
                            if otu_path:
                                try:
                                    fasta_text = otu_path.read_text(encoding="utf-8", errors="replace")
                                    combined_handle.write(prefix_fasta_ids(fasta_text, sample_id))
                                except Exception:
                                    pass
                            entry = _log_entry(utc_now_iso(), "OK", f"VSEARCH {sample_id}: OK")
                            logs.append(entry)
                            append_log(progress_path, entry)
                        else:
                            entry = _log_entry(utc_now_iso(), "ERR", f"VSEARCH {sample_id}: {vsearch_res.error_message or 'unknown error'}")
                            logs.append(entry)
                            append_log(progress_path, entry)
                            success = False

                combined_uc = run_config.output_dir / "analysis" / "combined_otus.uc"
                combined_global = run_config.output_dir / "analysis" / "combined_otus_global.fasta"
                sample_counts = {sample_record["internal_id"]: 1 for sample_record in next_samples}

                try:
                    subprocess.run(
                        [
                            VSEARCH_PATH,
                            "--cluster_size",
                            str(combined_fasta),
                            "--id",
                            str(run_config.vsearch_params.get("id", 0.97)),
                            "--centroids",
                            str(combined_global),
                            "--uc",
                            str(combined_uc),
                            "--threads",
                            str(run_config.vsearch_params.get("threads", 8)),
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                    )


                    matrix = _parse_vsearch_uc_to_matrix(combined_uc, sample_counts)
                    if matrix is not None:
                        matrix_path = run_config.output_dir / "analysis" / "otu_matrix.csv"
                        matrix.to_csv(str(matrix_path))
                        analysis_outputs["otu_matrix"] = str(matrix_path)
                        analysis_outputs["combined_otu_uc"] = str(combined_uc)

                except subprocess.CalledProcessError as exc:
                    success = False
                    entry = _log_entry(utc_now_iso(), "ERR", f"Global OTU clustering failed: {exc}")
                    logs.append(entry)
                    append_log(progress_path, entry)
                
                combined_fasta_seq_count = 0
                if combined_fasta.exists():
                    try:
                        combined_fasta_seq_count = sum(1 for line in combined_fasta.open(encoding="utf-8", errors="replace") if line.startswith(">"))
                    except Exception:
                        combined_fasta_seq_count = 0

                if combined_fasta_seq_count == 0:
                    success = False
                    entry = _log_entry(utc_now_iso(), "ERR", "Analysis: combined_fasta gerado com 0 sequências. Verifica se VSEARCH produziu OTUs.")
                    logs.append(entry)
                    append_log(progress_path, entry)
                    analysis_outputs["combined_fasta_vsearch"] = ""
                    analysis_outputs["combined_fasta"] = ""
                else:
                    analysis_outputs["combined_fasta_vsearch"] = str(combined_fasta)
                    analysis_outputs.setdefault("combined_fasta", str(combined_fasta))

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
                    success = False
                    entry = _log_entry(utc_now_iso(), "ERR", f"DADA2 failed: {dada2_result.error_message or 'unknown error'}")
                    logs.append(entry)
                    append_log(progress_path, entry)

            if run_vsearch_branch:
                vsearch_fasta = analysis_outputs.get("combined_fasta_vsearch") or None
                if vsearch_fasta:
                    vsearch_path = Path(vsearch_fasta)
                    vsearch_fasta = str(vsearch_path) if vsearch_path.exists() else None
                analysis_outputs["combined_fasta_vsearch"] = vsearch_fasta or ""
                analysis_outputs["combined_fasta"] = vsearch_fasta or ""
            elif run_dada2_branch:
                dada2_fasta = analysis_outputs.get("combined_fasta_dada2") or ""
                analysis_outputs["combined_fasta"] = dada2_fasta or ""
            else:
                analysis_outputs["combined_fasta"] = ""
                success = False
                entry = _log_entry(utc_now_iso(), "ERR", f"Analysis method inválido: {analysis_method}")
                logs.append(entry)
                append_log(progress_path, entry)

            analysis_record = {
                "signature": run_config.run_id,
                "status": "success" if success else "failure",
                "updated_at": utc_now_iso(),
                "outputs": analysis_outputs,
                "log_lines": logs.copy(),
            }
            progress_state["analysis_record"] = analysis_record
            _write_progress_state()
            stages["Analysis"] = analysis_record["status"]
            active_stage_results.append(analysis_record["status"] == "success")
            _publish()
        else:
            stages["Analysis"] = "skipped"
            progress_state["analysis_record"] = analysis_record
            _write_progress_state()
            _publish()

        # Databases
        databases_outputs = {}
        combined_fasta_value = analysis_outputs.get("combined_fasta")
        db_record = _skipped_record()
        if run_databases:
            stages["Databases"] = "running"
            _publish()
            if not combined_fasta_value:
                entry = _log_entry(utc_now_iso(), "ERR", "Databases: combined_fasta ausente")
                logs.append(entry)
                append_log(progress_path, entry)
                db_record = {
                    "signature": run_config.run_id,
                    "status": "failure",
                    "updated_at": utc_now_iso(),
                    "outputs": databases_outputs,
                    "log_lines": logs.copy(),
                }
                stages["Databases"] = "failure"
                active_stage_results.append(False)
                _publish()
            else:
                fasta_path = Path(combined_fasta_value)
                if not fasta_path.exists():
                    entry = _log_entry(utc_now_iso(), "ERR", f"Databases: combined_fasta inexistente: {fasta_path}")
                    logs.append(entry)
                    append_log(progress_path, entry)
                    db_record = {
                        "signature": run_config.run_id,
                        "status": "failure",
                        "updated_at": utc_now_iso(),
                        "outputs": databases_outputs,
                        "log_lines": logs.copy(),
                    }
                    stages["Databases"] = "failure"
                    active_stage_results.append(False)
                    _publish()
                elif analysis_record.get("status") not in ("success", "skipped"):
                    entry = _log_entry(utc_now_iso(), "ERR", "Databases: análise anterior falhou")
                    logs.append(entry)
                    append_log(progress_path, entry)
                    db_record = {
                        "signature": run_config.run_id,
                        "status": "failure",
                        "updated_at": utc_now_iso(),
                        "outputs": databases_outputs,
                        "log_lines": logs.copy(),
                    }
                    stages["Databases"] = "failure"
                    active_stage_results.append(False)
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
                    silva_config = SILVAConfig(db_path=db_paths.silva) if "silva" in db_sources and db_paths.silva and db_paths.silva.exists() else None
                    pr2_config = PR2Config(db_path=db_paths.pr2) if "pr2" in db_sources and db_paths.pr2 and db_paths.pr2.exists() else None
                    ncbi_config = (
                        NCBIConfig(
                            local_db_path=db_paths.blast if db_paths.blast and db_paths.blast.exists() else None,
                            use_remote=not (db_paths.blast and db_paths.blast.exists()),
                            evalue=float(run_config.vsearch_params.get("ncbi_evalue", 1e-10)),
                            max_target_seqs=int(run_config.vsearch_params.get("ncbi_maxhits", 5)),
                        )
                        if "ncbi" in db_sources
                        else None
                    )
                    bold_config = (
                        BOLDConfig(
                            db_path=db_paths.bold_dir,
                            db_source=int(run_config.bold_params.get("db", 1)),
                            mode=int(run_config.bold_params.get("mode", 1)),
                        )
                        if "bold" in db_sources and db_paths.bold_dir and db_paths.bold_dir.exists()
                        else None
                    )
                    tax_result = run_taxonomy_pipeline(
                        fasta_path,
                        run_config.output_dir / "databases",
                        silva_config=silva_config,
                        pr2_config=pr2_config,
                        ncbi_config=ncbi_config,
                        bold_config=bold_config,
                        custom_db_path=custom_db_path,
                        consensus_method=str(run_config.vsearch_params.get("db_consensus") or "lca"),
                        sintax_cutoff=sintax_cutoff,
                        threads=int(run_config.vsearch_params.get("threads", 8)),
                    )
                    databases_outputs["taxonomy"] = serialize_step_output(tax_result)
                    entry = _log_entry(
                        utc_now_iso(),
                        "OK" if tax_result.success else "ERR",
                        f"Taxonomy: {tax_result.error_message or 'OK'}" if not tax_result.success else "Taxonomy completed",
                    )
                    logs.append(entry)
                    append_log(progress_path, entry)
                    db_record = {
                        "signature": run_config.run_id,
                        "status": "success" if tax_result.success else "failure",
                        "updated_at": utc_now_iso(),
                        "outputs": databases_outputs,
                        "log_lines": logs.copy(),
                    }
                    progress_state["db_record"] = db_record
                    _write_progress_state()
                    stages["Databases"] = db_record["status"]
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
            fasta_for_phylo = (
                (db_record.get("outputs") or {}).get("taxonomy", {}).get("files", {}).get("annotated_fasta")
                or analysis_outputs.get("combined_fasta")
                or analysis_outputs.get("combined_fasta_vsearch")
                or analysis_outputs.get("combined_fasta_dada2")
            ) or None
            if fasta_for_phylo == "":
                fasta_for_phylo = None
            if not fasta_for_phylo:
                entry = _log_entry(utc_now_iso(), "ERR", "Phylo: input FASTA ausente")
                logs.append(entry)
                append_log(progress_path, entry)
                phylo_record = {
                    "signature": run_config.run_id,
                    "status": "failure",
                    "updated_at": utc_now_iso(),
                    "outputs": phylo_outputs,
                    "log_lines": logs.copy(),
                }
                stages["Phylo"] = "failure"
                active_stage_results.append(False)
                _publish()
            else:
                phylo_path = Path(fasta_for_phylo)
                if not phylo_path.exists():
                    entry = _log_entry(utc_now_iso(), "ERR", f"Phylo: input FASTA inexistente: {phylo_path}")
                    logs.append(entry)
                    append_log(progress_path, entry)
                    phylo_record = {
                        "signature": run_config.run_id,
                        "status": "failure",
                        "updated_at": utc_now_iso(),
                        "outputs": phylo_outputs,
                        "log_lines": logs.copy(),
                    }
                    stages["Phylo"] = "failure"
                    active_stage_results.append(False)
                    _publish()
                else:
                    try:
                        seq_count = sum(1 for _ in SeqIO.parse(str(phylo_path), "fasta"))
                    except Exception:
                        seq_count = 0
                    if seq_count < 2:
                        entry = _log_entry(utc_now_iso(), "ERR", f"Phylo: FASTA tem {seq_count} sequências (mínimo 2)")
                        logs.append(entry)
                        append_log(progress_path, entry)
                        phylo_record = {
                            "signature": run_config.run_id,
                            "status": "failure",
                            "updated_at": utc_now_iso(),
                            "outputs": phylo_outputs,
                            "log_lines": logs.copy(),
                        }
                        stages["Phylo"] = "failure"
                        active_stage_results.append(False)
                        _publish()
                        progress_state["phylo_record"] = phylo_record
                        _write_progress_state()

                    else:
                        entry = _log_entry(utc_now_iso(), "RUN", "Starting Phylo stage")
                        logs.append(entry)
                        append_log(progress_path, entry)
                        phylo_result = run_phylo_step(phylo_path, run_config, run_config.output_dir / "phylo")
                        phylo_outputs["phylo"] = {
                            "files": {k: str(v) for k, v in getattr(phylo_result, "files", {}).items()},
                            "metrics": getattr(phylo_result, "metrics", {}),
                        }
                        entry = _log_entry(
                            utc_now_iso(),
                            "OK" if phylo_result.success else "ERR",
                            f"Phylo: {phylo_result.error_message or 'OK'}" if not phylo_result.success else "Phylo completed",
                        )
                        logs.append(entry)
                        append_log(progress_path, entry)
                        phylo_record = {
                            "signature": run_config.run_id,
                            "status": "success" if phylo_result.success else "failure",
                            "updated_at": utc_now_iso(),
                            "outputs": phylo_outputs,
                            "log_lines": logs.copy(),
                        }
                        progress_state["phylo_record"] = phylo_record
                        _write_progress_state()
                        stages["Phylo"] = phylo_record["status"]
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
                "log_lines": logs.copy(),
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
                "log_lines": logs.copy(),
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


def register_callbacks(app: Any) -> None:
    @app.callback(
        Output("progress-interval", "disabled", allow_duplicate=True),
        Output("progress-overall-status", "children", allow_duplicate=True),
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
        progress_path = run_config.output_dir / "progress.json"
        thread = threading.Thread(target=_run_pipeline_thread, args=(run_config, run_config_data, progress_path), daemon=True)
        with _progress_lock:
            _pipeline_threads[run_config.run_id] = thread
        thread.start()
        return False, "Running..."

    @app.callback(
        Output("progress-stage-rows", "children"),
        Output("progress-logs", "children"),
        Output("progress-overall-status", "children", allow_duplicate=True),
        Output("progress-interval", "disabled", allow_duplicate=True),
        Output("store-qc-outputs", "data", allow_duplicate=True),
        Output("store-preprocess-outputs", "data", allow_duplicate=True),
        Output("store-analysis-outputs", "data", allow_duplicate=True),
        Output("store-databases-outputs", "data", allow_duplicate=True),
        Output("store-phylo-outputs", "data", allow_duplicate=True),
        Input("progress-interval", "n_intervals"),
        State("store-run-config", "data"),
        State("store-project", "data"),
        State("store-samples", "data"),
        prevent_initial_call=True,
    )
    def _poll_progress(n_intervals, run_config_data, project_data, samples_data):
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

        if done:
            return (
                stage_rows,
                log_block,
                overall_text,
                True,
                state.get("qc_record"),
                state.get("preprocess_record"),
                state.get("analysis_record"),
                state.get("db_record"),
                state.get("phylo_record"),
            )

        return stage_rows, log_block, overall_text, False, no_update, no_update, no_update, no_update, no_update