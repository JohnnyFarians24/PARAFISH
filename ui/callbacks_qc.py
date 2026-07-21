"""Callbacks related to the quality control stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, dash_table, dcc, html, no_update
from dash.exceptions import PreventUpdate

from core.models import AnalysisMode
from qc.fastqc_runner import run_fastqc
from ui.components.table import STYLE_CELL, data_table
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import latest_run_dir, resolve_project_data, run_config_from_store, serialize_step_output


def build_qc_summary(*args, **kwargs):
    """Build a small QC readiness summary from shared store data."""
    samples_data = kwargs.get("samples_data")
    step_outputs = kwargs.get("step_outputs") or {}
    project_data = resolve_project_data(kwargs.get("project_data") or {})
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


def _qc_metric_paths(qc_data, run_config_data, project_data=None) -> dict[str, str]:
    """Locate the QC tables the runner persisted.

    Prefers the paths recorded on the stage record; falls back to the conventional
    location under the run's output directory so runs recorded before the metrics
    layer existed still render.

    ``run_config_data`` is a per-session store: a fresh tab or a restarted app has it
    empty, so ``output_dir`` is gone even though the run is on disk. Fall back to the
    project's newest run (``latest_run_dir``, ordered by run_id) instead of blanking.
    """
    record = _normalize_store(qc_data)
    recorded = record.get("metrics_files") or {}
    paths = {key: value for key, value in recorded.items() if isinstance(value, str) and Path(value).exists()}
    if paths:
        return paths
    output_dir = (run_config_data or {}).get("output_dir")
    if not output_dir:
        run_dir = latest_run_dir(project_data) if project_data else None
        if run_dir is None:
            return {}
        output_dir = str(run_dir)
    from qc.qc_metrics import build_qc_metrics, qc_metrics_paths

    qc_dir = Path(output_dir) / "qc"
    existing = qc_metrics_paths(qc_dir)
    if existing:
        return {key: str(path) for key, path in existing.items()}
    # Runs recorded before the metrics layer existed have the FastQC zips but no
    # tables. Build them once and leave them on disk; this is a cache fill, not a
    # per-render recompute.
    try:
        return {key: str(path) for key, path in build_qc_metrics(qc_dir).items()}
    except Exception:
        return {}


_MODULE_STATUS_SYMBOL = {"pass": "●", "warn": "▲", "fail": "✕"}
_MODULE_STATUS_COLOR = {"pass": "var(--green)", "warn": "var(--amber)", "fail": "var(--red)"}


def _module_status_panel(status_csv: str):
    """Traffic-light grid of every FastQC module, per sample x read."""
    import pandas as pd

    try:
        status = pd.read_csv(status_csv)
    except Exception:
        return html.Div()
    if status.empty:
        return html.Div()

    modules = list(dict.fromkeys(status["module"]))
    header = html.Tr([html.Th("Sample / read"), *[html.Th(module, style={"fontSize": "11px"}) for module in modules]])
    body = []
    for (sample, read), group in status.groupby(["sample", "read"], sort=True):
        by_module = dict(zip(group["module"], group["status"]))
        cells = [html.Td(f"{sample} {read}", className="mono-cell")]
        for module in modules:
            value = str(by_module.get(module, "")).lower()
            cells.append(
                html.Td(
                    _MODULE_STATUS_SYMBOL.get(value, "–"),
                    title=f"{module}: {value or 'n/a'}",
                    style={"color": _MODULE_STATUS_COLOR.get(value, "var(--muted)"), "textAlign": "center", "fontWeight": "700"},
                )
            )
        body.append(html.Tr(cells))
    return html.Div(
        [
            html.Div("FastQC modules — ● pass · ▲ warn · ✕ fail", className="mini-note"),
            html.Div(html.Table([html.Thead(header), html.Tbody(body)], className="parafish-table"), style={"overflowX": "auto"}),
        ]
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
        Output("qc-plot-hero", "children"),
        Output("qc-plots", "children"),
        Input("store-qc-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_qc_plots(qc_out, run_config_data, project_data):
        """Split the hero (read depth) from the secondary curves.

        Read depth is the one plot that decides whether the run is usable, so it stays
        on the page. The per-base quality curves and the length histogram are what you
        consult *after* the headline tells you to look.
        """
        import pandas as pd

        from ui.plots import qc_depth_box_figure, qc_length_distribution_figure, qc_per_base_quality_figure

        def _graph(figure):
            return dcc.Graph(figure=figure, config={"displaylogo": False, "responsive": True})

        paths = _qc_metric_paths(qc_out, run_config_data, project_data)
        empty = html.Div("Run QC to see per-sample quality plots.", className="empty-state-sub")
        if not paths:
            return empty, empty
        try:
            hero = qc_depth_box_figure(pd.read_csv(paths["qc_overview_csv"])) if "qc_overview_csv" in paths else None
            secondary = [
                qc_per_base_quality_figure(pd.read_csv(paths["qc_per_base_quality_csv"])) if "qc_per_base_quality_csv" in paths else None,
                qc_length_distribution_figure(pd.read_csv(paths["qc_length_distribution_csv"])) if "qc_length_distribution_csv" in paths else None,
            ]
        except Exception as exc:
            error = html.Div(f"QC plots unavailable: {exc}", className="empty-state-sub")
            return error, error

        hero_child = _graph(hero) if hero is not None else html.Div("No read-depth data yet.", className="empty-state-sub")
        figures = [figure for figure in secondary if figure is not None]
        detail_child = (
            html.Div([_graph(fig) for fig in figures], style={"display": "flex", "flexDirection": "column", "gap": "12px"})
            if figures
            else html.Div("No FastQC metrics available yet.", className="empty-state-sub")
        )
        return hero_child, detail_child

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

        try:
            metrics_files = {key: str(path) for key, path in build_qc_metrics(run_config.output_dir / "qc").items()}
        except Exception as exc:
            metrics_files = {}
            log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"Falha ao construir métricas de QC: {exc}"})

        record = {
            "signature": signature,
            "status": "success" if all_success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": outputs,
            "metrics_files": metrics_files,
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
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_qc_sample_table(qc_data, run_config_data, project_data):
        import pandas as pd

        paths = _qc_metric_paths(qc_data, run_config_data, project_data)
        if "qc_overview_csv" not in paths:
            return html.Div("No QC results available yet.", className="empty-state-sub")
        try:
            overview = pd.read_csv(paths["qc_overview_csv"])
        except Exception as exc:
            return html.Div(f"QC overview unavailable: {exc}", className="empty-state-sub")
        if overview.empty:
            return html.Div("No QC results available yet.", className="empty-state-sub")

        columns = [
            {"name": "Sample", "id": "sample"},
            {"name": "Read", "id": "read"},
            {"name": "Total reads", "id": "total_reads", "type": "numeric", "format": {"specifier": ","}},
            {"name": "%GC", "id": "pct_gc"},
            {"name": "Mean length (bp)", "id": "mean_length"},
            {"name": "%duplication", "id": "pct_duplication"},
            {"name": "%adapter (max)", "id": "pct_adapter_max"},
            {"name": "%reads mean Q≥30", "id": "pct_reads_q30"},
            {"name": "FastQC fail", "id": "n_fail", "type": "numeric"},
            {"name": "FastQC warn", "id": "n_warn", "type": "numeric"},
        ]
        conditional = [
            {"if": {"filter_query": "{n_fail} > 0", "column_id": "n_fail"}, "color": "var(--red)", "fontWeight": "600"},
            {"if": {"filter_query": "{n_warn} > 0", "column_id": "n_warn"}, "color": "var(--amber)", "fontWeight": "600"},
            {"if": {"filter_query": "{pct_reads_q30} < 80", "column_id": "pct_reads_q30"}, "color": "var(--red)", "fontWeight": "600"},
            {"if": {"filter_query": "{pct_adapter_max} > 5", "column_id": "pct_adapter_max"}, "color": "var(--amber)"},
            {"if": {"filter_query": "{is_depth_outlier} = 1"}, "backgroundColor": "rgba(248,81,73,0.10)"},
        ]
        children = [
            data_table(
                overview.to_dict("records"),
                columns,
                conditional=conditional,
                sort_action="native",
                page_size=20,
                style_cell={**STYLE_CELL, "whiteSpace": "normal", "height": "auto"},
            )
        ]
        if int(overview["is_depth_outlier"].sum()) > 0:
            flagged = sorted(set(overview.loc[overview["is_depth_outlier"] == 1, "sample"]))
            children.append(
                html.Div(
                    f"Low-depth outlier(s), below half the median depth: {', '.join(flagged)}. "
                    "Consider removing them before alpha/beta diversity — they set the rarefaction floor.",
                    className="mini-note",
                    style={"color": "var(--red)"},
                )
            )
        if "qc_module_status_csv" in paths:
            children.append(_module_status_panel(paths["qc_module_status_csv"]))
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "10px"})
