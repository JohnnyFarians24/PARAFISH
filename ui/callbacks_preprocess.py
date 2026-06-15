"""Callbacks related to the preprocess stage."""

from __future__ import annotations

from typing import Any

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from preprocess.pear_runner import run_pear
from preprocess.trimmomatic_runner import run_trimmomatic
from ui.components.log_viewer import log_viewer
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import run_config_from_store, samples_from_store, serialize_step_output

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


def register_callbacks(app: Any) -> None:
    """Register preprocess callbacks in the Dash application."""
    @app.callback(
        Output("preprocess-summary", "children"),
        Output("preprocess-log-viewer", "children"),
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
            return _run_preprocess_auto(samples_data=samples_data, step_outputs={"qc": qc_record, "preprocess": existing}), log_viewer(existing.get("log_lines", [])), no_update

        if not samples:
            return _run_preprocess_auto(samples_data=samples_data, step_outputs={"qc": qc_record}), log_viewer([]), no_update

        if not stage_ready(qc_record):
            record = {
                "signature": signature,
                "status": "blocked" if stage_status(qc_record) == "failure" else "pending",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "next_samples": [],
                "log_lines": [
                    {"ts": utc_now_iso(), "level": "ERR" if stage_status(qc_record) == "failure" else "INF", "msg": "Preprocess aguardando QC concluído"}
                ],
                "blocked_by": "qc" if stage_status(qc_record) == "failure" else None,
            }
            return _run_preprocess_auto(samples_data=samples_data, step_outputs={"qc": qc_record, "preprocess": record}), log_viewer(record["log_lines"]), record

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
            return _run_preprocess_auto(samples_data=samples_data, step_outputs={"qc": qc_record, "preprocess": record}), log_viewer(record["log_lines"]), record

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

        record = {
            "signature": signature,
            "status": "success" if all_success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": outputs,
            "next_samples": next_samples,
            "log_lines": log_lines,
        }
        return _run_preprocess_auto(samples_data=samples_data, step_outputs={"qc": qc_record, "preprocess": record}), log_viewer(record["log_lines"]), record
