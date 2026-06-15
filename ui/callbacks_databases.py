"""Callbacks related to the databases stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from core.config import get_db_paths
from databases.bold_db import BOLDConfig
from databases.ncbi_db import NCBIConfig
from databases.pr2_db import PR2Config
from databases.silva_db import SILVAConfig
from databases.taxonomy_runner import run_taxonomy_pipeline
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import decode_upload_content, resolve_project_data, run_config_from_store, serialize_step_output


def _is_supported_custom_db_filename(filename: str) -> bool:
    name = (filename or "").lower()
    return name.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz"))


def run_db_classification(*args, **kwargs):
    """Build a compact database summary from shared store data."""
    project_data = kwargs.get("project_data") or {}
    samples_data = kwargs.get("samples_data") or {}
    step_outputs = kwargs.get("step_outputs") or {}
    stage_record = get_stage_record(step_outputs, "databases")
    sample_count = len(samples_data.get("samples", []))
    project_name = project_data.get("project_name") or "No project selected"
    status = stage_status(stage_record)
    return html.Div(
        [
            html.Div("Database routing", className="status-pill"),
            html.Div(f"Stage status: {status}", className="status-pill"),
            html.Div(f"Project: {project_name}", className="mini-note"),
            html.Div(f"{sample_count} sample group(s) awaiting taxonomy consensus", className="mini-note"),
        ],
        className="chip-list",
    )


def register_callbacks(app: Any) -> None:
    """Register database callbacks in the Dash application."""
    @app.callback(
        Output("custom-db-upload-status", "children"),
        Output("store-run-config", "data", allow_duplicate=True),
        Input("upload-custom-db", "contents"),
        State("upload-custom-db", "filename"),
        State("store-project", "data"),
        State("store-run-config", "data"),
        prevent_initial_call=True,
    )
    def _handle_custom_db_upload(contents, filename, project_data, run_config_data):
        if not contents or not filename:
            raise PreventUpdate
        if not _is_supported_custom_db_filename(filename):
            return html.Div(f"Unsupported format: {filename}. Use .fasta, .fa or .fasta.gz", style={"color": "red"}), no_update
        project = resolve_project_data(project_data)
        project_name = project.get("project_name") or "project"
        db_dir = Path("projects") / project_name / "custom_dbs"
        db_dir.mkdir(parents=True, exist_ok=True)
        dest = db_dir / Path(filename).name
        dest.write_bytes(decode_upload_content(contents))
        raw = dict(run_config_data or {})
        raw.setdefault("vsearch_params", {})["custom_db_path"] = str(dest)
        return html.Div(f"Uploaded: {filename} → {dest}", style={"color": "green"}), raw

    @app.callback(
        Output("databases-summary", "children"),
        Output("store-databases-outputs", "data"),
        Input("btn-run-databases", "n_clicks"),
        State("store-project", "data"),
        State("store-run-config", "data"),
        State("store-samples", "data"),
        State("store-analysis-outputs", "data"),
        State("store-databases-outputs", "data"),
        prevent_initial_call=True,
    )
    def _run_databases(n_clicks, project_data, run_config_data, samples_data, analysis_outputs, databases_outputs):
        if not n_clicks:
            raise PreventUpdate
        run_config = run_config_from_store(project_data, samples_data, run_config_data)
        analysis_record = _normalize_store(analysis_outputs)
        existing = _normalize_store(databases_outputs)
        signature = "|".join([run_config.run_id, run_config.branch.value, analysis_record.get("signature", "")])

        if existing.get("signature") == signature and stage_status(existing) in {"success", "failure", "blocked"}:
            return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": existing}), no_update

        if not stage_ready(analysis_record):
            record = {
                "signature": signature,
                "status": "blocked" if stage_status(analysis_record) == "failure" else "pending",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR" if stage_status(analysis_record) == "failure" else "INF", "msg": "Taxonomy aguardando análise concluída"}],
                "blocked_by": "analysis" if stage_status(analysis_record) == "failure" else None,
            }
            return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": record}), record

        combined_fasta_value = (analysis_record.get("outputs") or {}).get("combined_fasta")
        if not combined_fasta_value:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": "Análise não produziu FASTA combinado"}],
                "error_message": "Análise não produziu FASTA combinado",
            }
            return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": record}), record

        fasta_path = Path(combined_fasta_value)
        if not fasta_path.exists():
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": f"FASTA combinado nao encontrado: {fasta_path}"}],
                "error_message": f"FASTA combinado nao encontrado: {fasta_path}",
            }
            return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": record}), record
        db_paths = get_db_paths()
        custom_db_raw = (run_config.vsearch_params or {}).get("custom_db_path")
        custom_db_path = Path(custom_db_raw) if custom_db_raw and Path(custom_db_raw).exists() else None
        sintax_cutoff = float(run_config.vsearch_params.get("sintax_cutoff", 0.8))
        db_sources = set(
            (run_config_data or {}).get("db_sources")
            or (run_config_data or {}).get("vsearch_params", {}).get("db_sources")
            or []
        )
       
        silva_config = SILVAConfig(db_path=db_paths.silva) if "silva" in db_sources and db_paths.silva and db_paths.silva.exists() else None
        pr2_config = PR2Config(db_path=db_paths.pr2) if "pr2" in db_sources and db_paths.pr2 and db_paths.pr2.exists() else None
        bold_config = (
            BOLDConfig(
                db_path=db_paths.bold_dir,
                db_source=int(run_config.bold_params.get("db", 1)),
                mode=int(run_config.bold_params.get("mode", 1)),
            )
            if "bold" in db_sources and db_paths.bold_dir and db_paths.bold_dir.exists()
            else None
        )
        ncbi_config = (
            NCBIConfig(
                local_db_path=db_paths.blast,
                use_remote=not (db_paths.blast and db_paths.blast.exists()),
                evalue=float(run_config.vsearch_params.get("ncbi_evalue", 1e-10)),
                max_target_seqs=int(run_config.vsearch_params.get("ncbi_maxhits", 5)),
            )
            if "ncbi" in db_sources
            else None
        )

        log_lines = [{"ts": utc_now_iso(), "level": "RUN", "msg": "Starting taxonomy pipeline"}]
        taxonomy_result = run_taxonomy_pipeline(
			fasta_path,
			run_config.output_dir / "databases",
			silva_config=silva_config,
			pr2_config=pr2_config,
			ncbi_config=ncbi_config,
			bold_config=bold_config,
            custom_db_path=custom_db_path,
            sintax_cutoff=sintax_cutoff,
            consensus_method=str(run_config.vsearch_params.get("db_consensus") or "lca"),
            threads=int(run_config.vsearch_params.get("threads", 8)),
		)
        log_lines.append({"ts": utc_now_iso(), "level": "OK" if taxonomy_result.success else "ERR", "msg": taxonomy_result.error_message or "Taxonomy pipeline completed"})
        record = {
            "signature": signature,
            "status": "success" if taxonomy_result.success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": {"taxonomy": serialize_step_output(taxonomy_result)},
            "log_lines": log_lines,
            "error_message": taxonomy_result.error_message,
        }
        return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": record}), record
