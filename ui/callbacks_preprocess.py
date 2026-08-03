"""Callbacks related to the preprocess stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, dash_table, dcc, html, no_update
from dash.exceptions import PreventUpdate

from preprocess.pear_runner import run_pear
from preprocess.preprocess_metrics import build_preprocess_metrics
from preprocess.trimmomatic_runner import run_trimmomatic
from ui.components.log_viewer import log_viewer
from ui.components.table import data_table
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.plots import MIN_RETENTION_PCT
from ui.runtime_state import latest_run_dir, run_config_from_store, samples_from_store, serialize_step_output


def _preprocess_metric_paths(preprocess_data, run_config_data, project_data=None) -> dict[str, str]:
    """Locate the preprocess tables the runner persisted.

    Prefers the paths on the stage record; otherwise looks under the run's output
    directory, and as a last resort builds them once from the Trimmomatic/PEAR logs
    so runs recorded before this layer existed still render.

    ``run_config_data`` is a per-session store: a fresh tab or a restarted app has it
    empty, so ``output_dir`` is gone even though the run is still on disk. When that
    happens fall back to the project's newest run on disk (``latest_run_dir``) instead
    of returning ``{}`` and blanking the page — this is the disk path the Results page
    already takes, applied here.
    """
    record = _normalize_store(preprocess_data)
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
    from preprocess.preprocess_metrics import preprocess_metrics_paths

    pp_dir = Path(output_dir) / "preprocess"
    existing = preprocess_metrics_paths(pp_dir)
    if existing:
        return {key: str(path) for key, path in existing.items()}
    try:
        return {key: str(path) for key, path in build_preprocess_metrics(pp_dir).items()}
    except Exception:
        return {}


def _run_preprocess_auto(*args, **kwargs):
    """Build a compact preprocess summary from shared store data."""
    samples_data = kwargs.get("samples_data")
    step_outputs = kwargs.get("step_outputs") or {}
    stage_record = get_stage_record(step_outputs, "preprocess")
    sample_count = len((samples_data or {}).get("samples", []))
    status = stage_status(stage_record)
    return html.Div(
        [
            html.Div("Preprocess readiness", className="status-pill"),
            html.Div(f"Stage status: {status}", className="status-pill"),
            html.Div(f"{sample_count} sample group(s) queued for trimming and merging", className="mini-note"),
        ],
        className="chip-list",
    )


def _param_row(label: str, value: str, hint: str = ""):
    return html.Div(
        [
            html.Span(label, className="param-key"),
            html.Span(str(value), className="param-val"),
            html.Span(hint, className="param-hint") if hint else None,
        ],
        className="param-row",
    )


def register_callbacks(app: Any) -> None:

    @app.callback(
        Output("preprocess-params", "children"),
        Input("store-run-config", "data"),
    )
    def _render_preprocess_params(run_config_data):
        raw = run_config_data or {}
        trim = raw.get("trimmomatic_params") or {}
        pear = raw.get("pear_params") or {}
        window = f"{trim.get('window_size', 4)}:{trim.get('window_quality', 20)}"
        return html.Div(
            [
                html.Div("Trimmomatic — read trimming", className="param-section-title"),
                _param_row("Leading / Trailing", f"{trim.get('leading', 3)} / {trim.get('trailing', 3)}", "cut low-quality ends"),
                _param_row("Sliding window", window, "window:quality"),
                _param_row("Min length", trim.get("minlen", 50), "bp"),
                html.Div("PEAR — paired-end merging", className="param-section-title", style={"marginTop": "10px"}),
                _param_row("Min overlap", pear.get("min_overlap", 10), "bp"),
                _param_row("Min length", pear.get("min_length", 200), "bp"),
                _param_row("P-value", pear.get("p_value", 0.01), "assembly threshold"),
            ],
            className="param-grid",
        )

    @app.callback(
        Output("preprocess-log-viewer", "children"),
        Output("preprocess-summary", "children"),
        Input("store-preprocess-outputs", "data"),
        State("store-samples", "data"),
        State("store-qc-outputs", "data"),
    )
    def _sync_preprocess_from_store(preprocess_out, samples_data, qc_outputs):
        # The stage log only ever appeared when the user clicked "Run Preprocess"
        # on this page. A pipeline launched from Setup fills the store instead, so
        # the viewer sat on "No log lines yet" for the whole run.
        record = _normalize_store(preprocess_out)
        summary = _run_preprocess_auto(
            samples_data=samples_data,
            step_outputs={"qc": _normalize_store(qc_outputs), "preprocess": record},
        )
        return log_viewer(record.get("log_lines") or []), summary

    @app.callback(
        Output("preprocess-tracking-table", "children"),
        Input("store-preprocess-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_read_tracking(preprocess_out, run_config_data, project_data):
        import pandas as pd

        from preprocess.preprocess_metrics import failing_steps

        paths = _preprocess_metric_paths(preprocess_out, run_config_data, project_data)
        if "read_tracking_csv" not in paths:
            return html.Div("Run Preprocess to see the read-tracking waterfall.", className="empty-state-sub")
        try:
            tracking = pd.read_csv(paths["read_tracking_csv"])
        except Exception as exc:
            return html.Div(f"Read tracking unavailable: {exc}", className="empty-state-sub")
        if tracking.empty:
            return html.Div("No read-tracking data yet.", className="empty-state-sub")

        # The full waterfall, in pipeline order. "Filtered" and "Mapped" come from the
        # Analysis stage and are blank until it has run — they are what makes the row
        # add up, so they are shown rather than left in the CSV only.
        columns = [
            {"name": "Sample", "id": "sample"},
            {"name": "Mode", "id": "mode"},
            {"name": "Input", "id": "input_reads", "type": "numeric", "format": {"specifier": ","}},
            {"name": "With primer", "id": "reads_with_primer", "type": "numeric", "format": {"specifier": ","}},
            {"name": "% primer", "id": "pct_primer"},
            {"name": "Trimmed", "id": "trimmed_reads", "type": "numeric", "format": {"specifier": ","}},
            {"name": "% kept (trim)", "id": "pct_trimmed"},
            {"name": "Merged", "id": "assembled_reads", "type": "numeric", "format": {"specifier": ","}},
            {"name": "% kept (merge)", "id": "pct_assembled"},
            {"name": "Filtered", "id": "filtered_reads", "type": "numeric", "format": {"specifier": ","}},
            {"name": "% kept (filter)", "id": "pct_filtered"},
            {"name": "Mapped to OTU", "id": "mapped_reads", "type": "numeric", "format": {"specifier": ","}},
            {"name": "Unmapped", "id": "unmapped_reads", "type": "numeric", "format": {"specifier": ","}},
            {"name": "% of input", "id": "pct_mapped"},
        ]
        threshold = MIN_RETENTION_PCT
        conditional = [
            {"if": {"filter_query": f"{{pct_trimmed}} < {threshold}", "column_id": "pct_trimmed"}, "color": "var(--red)", "fontWeight": "700"},
            {"if": {"filter_query": f"{{pct_assembled}} < {threshold}", "column_id": "pct_assembled"}, "color": "var(--red)", "fontWeight": "700"},
            {"if": {"filter_query": f"{{pct_filtered}} < {threshold}", "column_id": "pct_filtered"}, "color": "var(--red)", "fontWeight": "700"},
            {"if": {"filter_query": f"{{pct_mapped}} < {threshold}", "column_id": "pct_mapped"}, "color": "var(--red)", "fontWeight": "700"},
        ]
        children = [
            data_table(
                tracking.to_dict("records"),
                columns,
                conditional=conditional,
                sort_action="native",
                page_size=20,
            )
        ]
        warnings = [
            f"{row['sample']}: {', '.join(steps)}"
            for row, steps in ((row, failing_steps(row)) for row in tracking.to_dict("records"))
            if steps
        ]
        if warnings:
            children.append(
                html.Div(
                    "Step(s) losing most of the reads — " + "; ".join(warnings),
                    className="mini-note",
                    style={"color": "var(--red)"},
                )
            )
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "8px"})

    @app.callback(
        Output("preprocess-stat-input", "children"),
        Output("preprocess-stat-surviving", "children"),
        Output("preprocess-stat-assembly", "children"),
        Input("store-preprocess-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_preprocess_stats(preprocess_out, run_config_data, project_data):
        import pandas as pd

        def _card(value, label, hint="", color=None):
            children = [html.Div(str(value), className="stat-value", style={"color": color} if color else {})]
            children.append(html.Div(label, className="stat-label"))
            if hint:
                children.append(html.Div(hint, className="stat-hint"))
            return children

        paths = _preprocess_metric_paths(preprocess_out, run_config_data, project_data)
        blank = (_card("—", "Input reads"), _card("—", "Surviving reads"), _card("—", "Assembly rate"))
        if "read_tracking_csv" not in paths:
            return blank
        try:
            tracking = pd.read_csv(paths["read_tracking_csv"])
        except Exception:
            return blank
        if tracking.empty:
            return blank

        total_input = int(tracking["input_reads"].sum())
        total_surviving = int(tracking["surviving_reads"].sum())
        pct_surviving = 100.0 * total_surviving / total_input if total_input else 0.0
        rates = pd.to_numeric(tracking["pct_assembled"], errors="coerce").dropna()
        if rates.empty:
            assembly = _card("n/a", "Assembly rate", "single-end run")
        else:
            worst = rates.min()
            assembly = _card(
                f"{rates.mean():.1f}%",
                "Assembly rate (mean)",
                f"worst sample {worst:.1f}%",
                color="var(--red)" if worst < MIN_RETENTION_PCT else None,
            )
        return (
            _card(f"{total_input:,}", "Input reads", "read pairs for PE"),
            _card(
                f"{total_surviving:,}",
                "Surviving reads",
                f"{pct_surviving:.1f}% of input",
                color="var(--red)" if pct_surviving < MIN_RETENTION_PCT else None,
            ),
            assembly,
        )

    @app.callback(
        Output("preprocess-plots", "children"),
        Input("store-preprocess-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_preprocess_plots(preprocess_out, run_config_data, project_data):
        import pandas as pd

        from ui.plots import assembly_rate_figure, fragment_length_figure, read_fate_figure

        paths = _preprocess_metric_paths(preprocess_out, run_config_data, project_data)
        if "read_tracking_csv" not in paths:
            return html.Div("Run Preprocess to see per-sample plots.", className="empty-state-sub")
        try:
            tracking = pd.read_csv(paths["read_tracking_csv"])
            figures = [read_fate_figure(tracking), assembly_rate_figure(tracking)]
            if "assembled_length_csv" in paths:
                figures.append(fragment_length_figure(pd.read_csv(paths["assembled_length_csv"])))
        except Exception as exc:
            return html.Div(f"Preprocess plots unavailable: {exc}", className="empty-state-sub")
        figures = [figure for figure in figures if figure is not None]
        if not figures:
            return html.Div("No preprocess metrics available yet.", className="empty-state-sub")
        return html.Div(
            [dcc.Graph(figure=fig, config={"displaylogo": False, "responsive": True}) for fig in figures],
            style={"display": "flex", "flexDirection": "column", "gap": "12px"},
        )
    """Register preprocess callbacks in the Dash application."""
    # The summary and log viewer are rendered by _sync_preprocess_from_store above,
    # so this callback only has to publish the stage record.
    @app.callback(
        Output("store-preprocess-outputs", "data"),
        Input("btn-run-preprocess", "n_clicks"),
        State("store-project", "data"),
        State("store-run-config", "data"),
        State("store-samples", "data"),
        State("store-qc-outputs", "data"),
        State("store-preprocess-outputs", "data"),
        prevent_initial_call=True,
    )
    def _run_preprocess(n_clicks, project_data, run_config_data, samples_data, qc_outputs, preprocess_outputs):
        if not n_clicks:
            raise PreventUpdate
        run_config = run_config_from_store(project_data, samples_data, run_config_data)
        samples = samples_from_store(samples_data)
        qc_record = _normalize_store(qc_outputs)
        existing = _normalize_store(preprocess_outputs)
        signature = "|".join([run_config.run_id, run_config.branch.value, qc_record.get("signature", ""), *[sample.internal_id for sample in samples]])

        if existing.get("signature") == signature and stage_status(existing) in {"success", "failure", "blocked"}:
            return no_update

        if not samples:
            return no_update

        if not stage_ready(qc_record):
            record = {
                "signature": signature,
                "status": "blocked" if stage_status(qc_record) == "failure" else "pending",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "next_samples": [],
                "log_lines": [
                    {"ts": utc_now_iso(), "level": "ERR" if stage_status(qc_record) == "failure" else "INF", "msg": "Preprocess waiting for QC to complete"}
                ],
                "blocked_by": "qc" if stage_status(qc_record) == "failure" else None,
            }
            return record

        input_errors = []
        for sample in run_config.samples:
            try:
                sample.paths.validate(sample.mode)
            except AssertionError as exc:
                input_errors.append(str(exc))
        if input_errors:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "next_samples": [],
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": error} for error in input_errors],
                "error_message": "; ".join(input_errors),
            }
            return record

        log_lines = [{"ts": utc_now_iso(), "level": "RUN", "msg": "Starting Trimmomatic stage"}]
        outputs = {}
        next_samples = []
        all_success = True

        for sample in run_config.samples:
            trim_result = run_trimmomatic(sample, run_config)
            outputs[sample.internal_id] = {"trimmomatic": serialize_step_output(trim_result)}
            if trim_result.success:
                log_lines.append({"ts": utc_now_iso(), "level": "OK", "msg": f"Trimmomatic completed: {sample.internal_id}"})
                next_sample = {
                    "internal_id": sample.internal_id,
                    "display_label": sample.display_label,
                    "mode": sample.mode.value,
                    "marker": sample.marker.value,
                    "is_control": sample.is_control,
                    "group": sample.group,
                    "r1_path": str(trim_result.files.get("r1_trimmed", sample.paths.r1)),
                    "r2_path": str(trim_result.files.get("r2_trimmed")) if trim_result.files.get("r2_trimmed") else None,
                }
                next_samples.append(next_sample)
            else:
                all_success = False
                log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"Trimmomatic failed: {sample.internal_id} -> {trim_result.error_message or 'unknown error'}"})
                continue

            if sample.mode.value == "PE":
                pear_result = run_pear(sample, run_config, trim_result)
                outputs[sample.internal_id]["pear"] = serialize_step_output(pear_result)
                if pear_result.success:
                    log_lines.append({"ts": utc_now_iso(), "level": "OK", "msg": f"PEAR completed: {sample.internal_id}"})
                else:
                    all_success = False
                    log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"PEAR failed: {sample.internal_id} -> {pear_result.error_message or 'unknown error'}"})

        try:
            metrics_files = {key: str(path) for key, path in build_preprocess_metrics(run_config.output_dir / "preprocess").items()}
        except Exception as exc:
            metrics_files = {}
            log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"Failed to build preprocess metrics: {exc}"})

        record = {
            "signature": signature,
            "status": "success" if all_success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": outputs,
            "next_samples": next_samples,
            "metrics_files": metrics_files,
            "log_lines": log_lines,
        }
        return record
