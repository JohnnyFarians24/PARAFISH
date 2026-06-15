"""Callbacks related to the analysis stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from analysis.dada2_runner import run_dada2
from analysis.vsearch_runner import run_vsearch
from core.models import AnalysisBranch, AnalysisMode, Marker, Sample, SamplePaths
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import resolve_default_sintax_database, run_config_from_store, samples_from_store, serialize_step_output, build_sample_from_record


def run_analysis(*args, **kwargs):
    """Build a compact analysis summary from shared store data."""
    samples_data = kwargs.get("samples_data")
    step_outputs = kwargs.get("step_outputs") or {}
    stage_record = get_stage_record(step_outputs, "analysis")
    sample_count = len((samples_data or {}).get("samples", []))
    status = stage_status(stage_record)
    return html.Div(
        [
            html.Div("Analysis readiness", className="status-pill"),
            html.Div(f"Stage status: {status}", className="status-pill"),
            html.Div(f"{sample_count} sample group(s) available for VSEARCH/DADA2", className="mini-note"),
        ],
        className="chip-list",
    )





def register_callbacks(app: Any) -> None:

    @app.callback(
        Output("analysis-stat-otus", "children"),
        Output("analysis-stat-asvs", "children"),
        Output("analysis-stat-consensus", "children"),
        Input("store-analysis-outputs", "data"),
    )
    def _update_analysis_stats(analysis_data):
        from core.utils import _count_fasta_seqs
        data = _normalize_store(analysis_data)
        outputs = data.get("outputs") or {}

        def _card(value, label):
            return [html.Div(str(value), className="stat-value"), html.Div(label, className="stat-label")]

        combined_fasta = outputs.get("combined_fasta") or ""
        otu_count = "—"
        if combined_fasta and Path(combined_fasta).exists():
            otu_count = _count_fasta_seqs(Path(combined_fasta))

        dada2_out = outputs.get("dada2") or {}
        asv_fasta = (dada2_out.get("files") or {}).get("asv_sequences_fasta") or ""
        asv_count = "—"
        if asv_fasta and Path(asv_fasta).exists():
            asv_count = _count_fasta_seqs(Path(asv_fasta))

        tax_out = outputs.get("taxonomy") or {}
        consensus = (tax_out.get("files") or {}).get("taxonomy_consensus_csv") or ""
        consensus_label = "available" if consensus and Path(consensus).exists() else "—"

        return _card(otu_count, "OTUs"), _card(asv_count, "ASVs"), _card(consensus_label, "Consensus")

    @app.callback(
        Output("analysis-summary", "children"),
        Output("store-analysis-outputs", "data"),
        Input("btn-run-analysis", "n_clicks"),
        State("store-project", "data"),
        State("store-run-config", "data"),
        State("store-samples", "data"),
        State("store-preprocess-outputs", "data"),
        State("store-analysis-outputs", "data"),
        prevent_initial_call=True,
    )
    def _run_analysis(n_clicks, project_data, run_config_data, samples_data, preprocess_outputs, analysis_outputs):
        if not n_clicks:
            raise PreventUpdate
        run_config = run_config_from_store(project_data, samples_data, run_config_data)
        preprocess_record = preprocess_outputs or {}
        existing = analysis_outputs or {}
        signature = "|".join([run_config.run_id, run_config.branch.value, preprocess_record.get("signature", "")])

        if existing.get("signature") == signature and stage_status(existing) in {"success", "failure", "blocked"}:
            return run_analysis(samples_data=samples_data, step_outputs={"preprocess": preprocess_record, "analysis": existing}), no_update

        if not stage_ready(preprocess_record):
            record = {
                "signature": signature,
                "status": "blocked" if stage_status(preprocess_record) == "failure" else "pending",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR" if stage_status(preprocess_record) == "failure" else "INF", "msg": "Analysis aguardando preprocess concluído"}],
                "blocked_by": "preprocess" if stage_status(preprocess_record) == "failure" else None,
            }
            return run_analysis(samples_data=samples_data, step_outputs={"preprocess": preprocess_record, "analysis": record}), record

        input_errors = []
        if run_config.branch == AnalysisBranch.VSEARCH:
            for sample_record in (preprocess_record.get("next_samples") or []):
                sample = _build_sample_from_paths(sample_record)
                try:
                    sample.paths.validate(sample.mode)
                except AssertionError as exc:
                    input_errors.append(str(exc))
        else:
            for sample_record in (preprocess_record.get("next_samples") or []):
                sample = _build_sample_from_paths(sample_record)
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
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": error} for error in input_errors],
                "error_message": "; ".join(input_errors),
            }
            return run_analysis(samples_data=samples_data, step_outputs={"preprocess": preprocess_record, "analysis": record}), record

        log_lines = [{"ts": utc_now_iso(), "level": "RUN", "msg": f"Analysis branch: {run_config.branch.value}"}]
        outputs = {}
        success = True

        if run_config.branch == AnalysisBranch.VSEARCH:
            combined_fasta = run_config.output_dir / "analysis" / "combined_otus.fasta"
            combined_fasta.parent.mkdir(parents=True, exist_ok=True)
            with combined_fasta.open("w", encoding="utf-8") as combined_handle:
                for sample_record in (preprocess_record.get("next_samples") or []):
                    sample_id = sample_record["internal_id"]
                    sample_marker = str(sample_record.get("marker") or "")
                    pear_record = (preprocess_record.get("outputs") or {}).get(sample_id, {}).get("pear")
                    if not pear_record:
                        success = False
                        log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"Missing PEAR output for {sample_id}"})
                        continue
                    assembled_value = (pear_record.get("files") or {}).get("assembled")
                    if not assembled_value:
                        success = False
                        log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"PEAR assembled path missing for {sample_id}"})
                        continue
                    assembled_path = Path(assembled_value)
                    if not assembled_path.exists():
                        success = False
                        log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"PEAR assembled file missing for {sample_id}"})
                        continue
                    sample = _build_sample_from_paths(sample_record)
                    db_sintax = resolve_default_sintax_database(sample_marker)
                    vsearch_result = run_vsearch(assembled_path, run_config, run_config.output_dir / "analysis" / sample.internal_id, db_sintax=Path(db_sintax) if db_sintax else None)
                    outputs[sample.internal_id] = {"vsearch": serialize_step_output(vsearch_result)}
                    if vsearch_result.success:
                        otu_path = vsearch_result.files.get("otus_fasta")
                        if otu_path and otu_path.exists():
                            combined_handle.write(otu_path.read_text(encoding="utf-8", errors="replace"))
                        log_lines.append({"ts": utc_now_iso(), "level": "OK", "msg": f"VSEARCH completed: {sample.internal_id}"})
                    else:
                        success = False
                        log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"VSEARCH failed: {sample.internal_id} -> {vsearch_result.error_message or 'unknown error'}"})
            outputs["combined_fasta"] = str(combined_fasta)
        else:
            analysis_samples = [_build_sample_from_paths(sample_record) for sample_record in (preprocess_record.get("next_samples") or [])]
            dada2_result = run_dada2(analysis_samples, run_config, run_config.output_dir / "analysis" / "dada2")
            outputs["dada2"] = serialize_step_output(dada2_result)
            if dada2_result.success:
                outputs["combined_fasta"] = str(dada2_result.files.get("asv_sequences_fasta"))
                log_lines.append({"ts": utc_now_iso(), "level": "OK", "msg": "DADA2 completed"})
            else:
                success = False
                log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"DADA2 failed -> {dada2_result.error_message or 'unknown error'}"})

        record = {
            "signature": signature,
            "status": "success" if success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": outputs,
            "log_lines": log_lines,
        }
        return run_analysis(samples_data=samples_data, step_outputs={"preprocess": preprocess_record, "analysis": record}), record
