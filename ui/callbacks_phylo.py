"""Callbacks related to the phylogenetics stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio import SeqIO
from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from phylo.phylo_runner import run_phylo as run_phylo_step
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import run_config_from_store


def _build_phylo_summary(**kwargs):
    """Build a compact phylogenetics summary from shared store data."""
    step_outputs = kwargs.get("step_outputs") or {}
    sample_count = len((kwargs.get("samples_data") or {}).get("samples", []))
    stage_record = get_stage_record(step_outputs, "phylo")
    status = stage_status(stage_record)
    return html.Div(
        [
            html.Div("Phylogenetics readiness", className="status-pill"),
            html.Div(f"Stage status: {status}", className="status-pill"),
            html.Div(f"{sample_count} sample group(s) ready for tree inference", className="mini-note"),
            html.Div(f"{len(step_outputs)} stored output group(s) available for tree building", className="mini-note"),
        ],
        className="chip-list",
    )


def register_callbacks(app: Any) -> None:
    """Register phylogenetics callbacks in the Dash application."""
    @app.callback(
        Output("phylo-summary", "children"),
        Output("store-phylo-outputs", "data"),
        Input("btn-run-phylo", "n_clicks"),
        State("store-project", "data"),
        State("store-run-config", "data"),
        State("store-samples", "data"),
        State("store-analysis-outputs", "data"),
        State("store-databases-outputs", "data"),
        State("store-phylo-outputs", "data"),
        prevent_initial_call=True,
    )
    def _run_phylo(n_clicks, project_data, run_config_data, samples_data, analysis_outputs, databases_outputs, phylo_outputs):
        if not n_clicks:
            raise PreventUpdate
        run_config = run_config_from_store(project_data, samples_data, run_config_data)
        analysis_record = _normalize_store(analysis_outputs)
        databases_record = _normalize_store(databases_outputs)
        existing = _normalize_store(phylo_outputs)
        signature = "|".join([run_config.run_id, run_config.phylo_mode.value if run_config.phylo_mode else "", databases_record.get("signature", ""), analysis_record.get("signature", "")])

        if existing.get("signature") == signature and stage_status(existing) in {"success", "failure", "blocked"}:
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": existing}), no_update

        if stage_status(databases_record) in {"failure", "blocked"}:
            blocker = "databases"
            record = {
                "signature": signature,
                "status": "blocked" if blocker else "pending",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "INF", "msg": "Phylo aguardando entradas concluídas"}],
                "blocked_by": blocker,
            }
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

        if not stage_ready(databases_record):
            record = {
                "signature": signature,
                "status": "pending",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "INF", "msg": "Phylo aguardando taxonomia concluída"}],
            }
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

        fasta_value = None
        if stage_ready(databases_record):
            fasta_value = (databases_record.get("outputs") or {}).get("taxonomy", {}).get("files", {}).get("annotated_fasta")
        if not fasta_value:
            fasta_value = (analysis_record.get("outputs") or {}).get("combined_fasta")
        if not fasta_value:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": "Nenhum FASTA disponível para phylo"}],
                "error_message": "Nenhum FASTA disponível para phylo",
            }
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

        fasta_path = Path(fasta_value)
        if not fasta_path.exists():
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": f"FASTA nao encontrado: {fasta_path}"}],
                "error_message": f"FASTA nao encontrado: {fasta_path}",
            }
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

        if fasta_path.stat().st_size == 0:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": f"FASTA vazio: {fasta_path}"}],
                "error_message": f"FASTA vazio: {fasta_path}",
            }
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

        try:
            seq_count = sum(1 for _ in SeqIO.parse(str(fasta_path), "fasta"))
        except Exception as exc:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": f"Falha a ler FASTA para phylo: {exc}"}],
                "error_message": f"Falha a ler FASTA para phylo: {exc}",
            }
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

        if seq_count < 2:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": f"FASTA tem {seq_count} sequências; mínimo 2 para phylo"}],
                "error_message": f"FASTA tem {seq_count} sequências; mínimo 2 para phylo",
            }
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

        phylo_result = run_phylo_step(fasta_path, run_config, run_config.output_dir / "phylo")
        record = {
            "signature": signature,
            "status": "success" if phylo_result.success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": {"phylo": {"files": {key: str(value) for key, value in phylo_result.files.items()}, "metrics": phylo_result.metrics}},
            "log_lines": [{"ts": utc_now_iso(), "level": "OK" if phylo_result.success else "ERR", "msg": phylo_result.error_message or "Phylo pipeline completed"}],
            "error_message": phylo_result.error_message,
        }
        return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

    @app.callback(
        Output("phylo-pipeline-bar", "children"),
        Input("store-phylo-outputs", "data"),
    )
    def _render_phylo_viewer(phylo_data):
        from core.config import get_runs_dir
        data = _normalize_store(phylo_data)
        outputs = data.get("outputs") or {}
        phylo_out = outputs.get("phylo") or {}
        files = phylo_out.get("files") or {}
        runs_root = Path(get_runs_dir())

        items = []
        for key, label in [("best_tree", "Best Tree (.nwk)"), ("support_tree", "Support Tree"), ("bootstraps", "Bootstraps"), ("best_model", "Best Model")]:
            path_str = files.get(key)
            if not path_str:
                continue
            p = Path(path_str)
            if not p.exists():
                continue
            try:
                rel = p.relative_to(runs_root)
                href = f"/runs/{rel.as_posix()}"
            except ValueError:
                href = "#"
            items.append(html.Div(html.A(f"⬇ {label}", href=href, target="_blank", className="btn-parafish")))

        if not items:
            return html.Div("No phylogenetic outputs available yet.", className="empty-state-sub")
        return html.Div(items, style={"display": "flex", "flexDirection": "column", "gap": "8px"})