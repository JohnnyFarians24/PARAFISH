"""Callbacks related to the quality control stage."""

from __future__ import annotations

from typing import Any

from dash import Input, Output, State, dash_table, html, no_update
from dash.exceptions import PreventUpdate

from core.models import AnalysisMode
from qc.fastqc_runner import run_fastqc
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import run_config_from_store, serialize_step_output


def build_qc_summary(*args, **kwargs):
    """Build a small QC readiness summary from shared store data."""
    samples_data = kwargs.get("samples_data")
    step_outputs = kwargs.get("step_outputs") or {}
    project_data = kwargs.get("project_data") or {}
    stage_record = get_stage_record(step_outputs, "qc")
    sample_count = len((samples_data or {}).get("samples", []))
    output_count = len((stage_record or {}).get("outputs", {}))
    project_name = project_data.get("project_name") or "No project selected"
    status = stage_status(stage_record)
    return html.Div(
        [
            html.Div(f"Project: {project_name}", className="status-pill"),
            html.Div(f"QC status: {status}", className="status-pill"),
            html.Div(f"{sample_count} sample group(s) ready for QC", className="mini-note"),
            html.Div(f"{output_count} finished step output(s) cached", className="mini-note"),
        ],
        className="chip-list",
    )


def _qc_preflight(run_config) -> list[str]:
    errors = []
    seen_ids = set()
    for sample in run_config.samples:
        if sample.internal_id in seen_ids:
            errors.append(f"Sample ID duplicado: {sample.internal_id}")
        seen_ids.add(sample.internal_id)
        try:
            sample.paths.validate(sample.mode)
        except AssertionError as exc:
            errors.append(str(exc))
        if sample.mode == AnalysisMode.PE and sample.paths.r2 is None:
            errors.append(f"Modo PE requer R2: {sample.internal_id}")
    return errors


def register_callbacks(app: Any) -> None:
    """Register QC callbacks in the Dash application."""
    @app.callback(
        Output("qc-stat-total", "children"),
        Output("qc-stat-pass", "children"),
        Output("qc-stat-reads", "children"),
        Output("qc-stat-gc", "children"),
        Input("store-qc-outputs", "data"),
    )
    def update_qc_stats(qc_data):
        def _card(value, label, color=None):
            style = {"color": color} if color else {}
            return [html.Div(str(value), className="stat-value", style=style), html.Div(label, className="stat-label")]
        qc_data = _normalize_store(qc_data)

        status = str(qc_data.get("status") or "pending")
        outputs = _normalize_store(qc_data.get("outputs"))
        if status in {"failure", "blocked"}:
            return (
                _card("ERR", "Total samples", color="var(--red)"),
                _card("ERR", "Pass rate", color="var(--red)"),
                _card("ERR", "Mean reads", color="var(--red)"),
                _card("ERR", "Mean GC", color="var(--red)"),
            )

        if not qc_data:
            return (
                _card("—", "Total samples"),
                _card("—", "Pass rate"),
                _card("—", "Mean reads"),
                _card("—", "Mean GC"),
            )

        total = len(outputs)
        passed = sum(1 for output in outputs.values() if output.get("success"))
        pass_pct = f"{(100 * passed) // total}%" if total else "—"
        all_reads = [output.get("metrics", {}).get("total_sequences_r1", 0) for output in outputs.values()]
        mean_reads = f"{int(sum(all_reads) / len(all_reads)):,}" if all_reads else "—"
        all_gc = [output.get("metrics", {}).get("pct_gc_r1", 0) for output in outputs.values()]
        mean_gc = f"{(sum(all_gc) / len(all_gc)):.0f}%" if all_gc else "—"

        return (
            _card(total, "Total samples"),
            _card(pass_pct, "Pass rate", color="var(--green)" if passed == total else "var(--amber)"),
            _card(mean_reads, "Mean reads"),
            _card(mean_gc, "Mean GC"),
        )

    @app.callback(
        Output("qc-summary", "children"),
        Output("store-qc-outputs", "data"),
        Input("btn-run-qc", "n_clicks"),
        State("store-project", "data"),
        State("store-run-config", "data"),
        State("store-samples", "data"),
        State("store-qc-outputs", "data"),
        prevent_initial_call=True,
    )
    def _run_qc(n_clicks, project_data, run_config_data, samples_data, qc_outputs):
        if not n_clicks:
            raise PreventUpdate
        run_config = run_config_from_store(project_data, samples_data, run_config_data)
        if not run_config.samples:
            return build_qc_summary(project_data=project_data, samples_data=samples_data, step_outputs={"qc": {"status": "pending"}}), no_update
        signature = "|".join(
            [
                run_config.run_id,
                run_config.branch.value,
                *[
                    f"{sample.internal_id}:{sample.paths.r1}:{sample.paths.r2 or ''}:{sample.mode.value}"
                    for sample in run_config.samples
                ],
            ]
        )

        errors = _qc_preflight(run_config)
        if errors:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "errors": errors,
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": error} for error in errors],
            }
            return build_qc_summary(project_data=project_data, samples_data=samples_data, step_outputs={"qc": record}), record

        outputs = {}
        log_lines = []
        all_success = True
        for sample in run_config.samples:
            log_lines.append({"ts": utc_now_iso(), "level": "RUN", "msg": f"FastQC start: {sample.internal_id}"})
            result = run_fastqc(sample, run_config.output_dir / "qc" / sample.internal_id, threads=int(run_config.vsearch_params.get("threads", 4)))
            outputs[sample.internal_id] = serialize_step_output(result)
            if result.success:
                log_lines.append({"ts": utc_now_iso(), "level": "OK", "msg": f"FastQC completed: {sample.internal_id}"})
            else:
                all_success = False
                log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"FastQC failed: {sample.internal_id} -> {result.error_message or 'unknown error'}"})

        record = {
            "signature": signature,
            "status": "success" if all_success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": outputs,
            "log_lines": log_lines,
        }
        return build_qc_summary(project_data=project_data, samples_data=samples_data, step_outputs={"qc": record}), record

    @app.callback(
        Output("qc-summary", "children", allow_duplicate=True),
        Input("store-qc-outputs", "data"),
        State("store-project", "data"),
        State("store-samples", "data"),
        prevent_initial_call=True,
    )
    def _update_qc_summary_from_store(qc_data, project_data, samples_data):
        qc_data = _normalize_store(qc_data)
        if not qc_data or qc_data.get("status") in (None, "pending", "skipped"):
            raise PreventUpdate
        return build_qc_summary(
            project_data=project_data,
            samples_data=samples_data,
            step_outputs={"qc": qc_data},
        )

    @app.callback(
        Output("qc-sample-table", "children"),
        Input("store-qc-outputs", "data"),
    )
    def _render_qc_sample_table(qc_data):
        qc_data = _normalize_store(qc_data)
        outputs = _normalize_store(qc_data.get("outputs")) if qc_data else {}
        if not outputs:
            return html.Div("No QC results available yet.", className="empty-state-sub")

        rows = []
        for sample_id, sample_data in outputs.items():
            metrics = sample_data.get("metrics") or {}
            rows.append(
                {
                    "Sample": sample_id,
                    "Total Reads R1": metrics.get("total_sequences_r1", "—"),
                    "Total Reads R2": metrics.get("total_sequences_r2", "—"),
                    "GC% R1": metrics.get("pct_gc_r1", "—"),
                    "GC% R2": metrics.get("pct_gc_r2", "—"),
                    "Failed Modules": ", ".join(metrics.get("failed_modules") or []),
                    "Status": "success" if sample_data.get("success") else "failure",
                }
            )

        table = dash_table.DataTable(
            data=rows,
            columns=[
                {"name": "Sample", "id": "Sample"},
                {"name": "Total Reads R1", "id": "Total Reads R1"},
                {"name": "Total Reads R2", "id": "Total Reads R2"},
                {"name": "GC% R1", "id": "GC% R1"},
                {"name": "GC% R2", "id": "GC% R2"},
                {"name": "Failed Modules", "id": "Failed Modules"},
                {"name": "Status", "id": "Status"},
            ],
            page_size=20,
            style_table={"overflowX": "auto"},
            style_cell={"padding": "8px", "whiteSpace": "normal", "height": "auto"},
            style_header={"fontWeight": "600"},
        )
        return table
