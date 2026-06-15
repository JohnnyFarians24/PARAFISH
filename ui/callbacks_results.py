"""Populate the Results page from stored step outputs and expose download links."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, html
from dash.exceptions import PreventUpdate

from core.config import get_runs_dir
from ui.pipeline_state import _normalize_store


def _collect_files_from_outputs(outputs: dict) -> list[str]:
    files = []

    def walk(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k in {"combined_fasta", "asv_sequences_fasta", "tree_nwk", "newick"} and isinstance(v, str):
                    files.append(v)
                elif k == "files" and isinstance(v, dict):
                    for item in v.values():
                        if isinstance(item, str):
                            files.append(item)
                else:
                    walk(v)

    walk(outputs or {})
    return files


def _safe_relative(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def register_callbacks(app: Any) -> None:

    @app.callback(
        Output("store-qc-outputs", "data", allow_duplicate=True),
        Output("store-preprocess-outputs", "data", allow_duplicate=True),
        Output("store-analysis-outputs", "data", allow_duplicate=True),
        Output("store-databases-outputs", "data", allow_duplicate=True),
        Output("store-phylo-outputs", "data", allow_duplicate=True),
        Output("store-run-config", "data", allow_duplicate=True),
        Output("url", "pathname", allow_duplicate=True),
        Input("btn-new-run", "n_clicks"),
        prevent_initial_call=True,
    )
    def _new_run(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        return {}, {}, {}, {}, {}, {}, "/setup"

    @app.callback(
        Output("results-stage-selector", "options"),
        Input("store-qc-outputs", "data"),
        Input("store-preprocess-outputs", "data"),
        Input("store-analysis-outputs", "data"),
        Input("store-databases-outputs", "data"),
        Input("store-phylo-outputs", "data"),
    )
    def _populate_stage_selector(qc, pre, analysis, db, phylo):
        stage_map = {"qc": qc, "preprocess": pre, "analysis": analysis, "databases": db, "phylo": phylo}
        return [
            {"label": name.upper(), "value": name}
            for name, data in stage_map.items()
            if (data or {}).get("status") == "success"
        ]

    @app.callback(
        Output("results-sample-selector", "options"),
        Input("results-stage-selector", "value"),
        State("store-qc-outputs", "data"),
        State("store-preprocess-outputs", "data"),
        State("store-analysis-outputs", "data"),
        State("store-databases-outputs", "data"),
        State("store-phylo-outputs", "data"),
    )
    def _populate_sample_selector(stage, qc, pre, analysis, db, phylo):
        if not stage:
            return []
        stage_map = {"qc": qc, "preprocess": pre, "analysis": analysis, "databases": db, "phylo": phylo}
        record = _normalize_store(stage_map.get(stage))
        outputs = record.get("outputs") or {}
        sample_options = [
            {"label": sid, "value": sid}
            for sid in outputs
            if sid not in ("combined_fasta", "dada2", "taxonomy")
        ]

        if stage in ("analysis", "databases", "phylo") and not sample_options:
            return [{"label": "Overview", "value": "__ALL__"}]
        
        return sample_options

    @app.callback(
        Output("results-viewer", "children"),
        Input("results-stage-selector", "value"),
        Input("results-sample-selector", "value"),
        State("store-qc-outputs", "data"),
        State("store-preprocess-outputs", "data"),
        State("store-analysis-outputs", "data"),
        State("store-databases-outputs", "data"),
        State("store-phylo-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_viewer(stage, sample_id, qc, pre, analysis, db, phylo, run_config_data):
        if sample_id in (None, "", "__ALL__"):
            sample_id = None
            
        if not stage:
            return html.Div("Select a stage to preview outputs.", className="empty-state-sub")
        stage_map = {"qc": qc, "preprocess": pre, "analysis": analysis, "databases": db, "phylo": phylo}
        record = _normalize_store(stage_map.get(stage))
        outputs = record.get("outputs") or {}
        runs_root = Path(get_runs_dir())

        if stage == "qc" and sample_id:
            sample_out = outputs.get(sample_id) or {}
            files = sample_out.get("files") or {}
            html_path = next((v for v in files.values() if str(v).endswith(".html")), None)
            if html_path:
                try:
                    rel = Path(html_path).relative_to(runs_root)
                    href = f"/runs/{rel.as_posix()}"
                    return html.Iframe(src=href, style={"width": "100%", "height": "600px", "border": "none"})
                except ValueError:
                    pass
            return html.Div("QC HTML not available.", className="empty-state-sub")

        if stage == "preprocess" and sample_id:
            sample_out = outputs.get(sample_id) or {}
            trim_files = (sample_out.get("trimmomatic") or {}).get("files") or {}
            pear_files = (sample_out.get("pear") or {}).get("files") or {}
            items = []
            for label, path in {**trim_files, **pear_files}.items():
                if path and Path(str(path)).exists():
                    try:
                        rel = Path(str(path)).relative_to(runs_root)
                        items.append(html.Div(html.A(f"{label}: {Path(path).name}", href=f"/runs/{rel.as_posix()}", target="_blank")))
                    except ValueError:
                        pass
            return html.Div(items) if items else html.Div("No preprocess outputs found.", className="empty-state-sub")

        if stage == "analysis":
            vsearch_out = outputs.get(sample_id, {}).get("vsearch") if sample_id else None
            otu_path = (vsearch_out.get("files") or {}).get("otus_fasta") if vsearch_out else outputs.get("combined_fasta")
            if otu_path and Path(str(otu_path)).exists():
                from core.utils import _count_fasta_seqs

                count = _count_fasta_seqs(Path(otu_path))
                try:
                    rel = Path(str(otu_path)).relative_to(runs_root)
                    href = f"/runs/{rel.as_posix()}"
                except ValueError:
                    href = "#"
                return html.Div([
                    html.Div(f"{count} OTUs/sequences", className="stat-value"),
                    html.A("Download FASTA", href=href, target="_blank"),
                ])
            return html.Div("No analysis FASTA available.", className="empty-state-sub")

        if stage == "databases":
            tax_record = outputs.get("taxonomy") or {}
            consensus_path = (tax_record.get("files") or {}).get("taxonomy_consensus_csv")
            if consensus_path and Path(consensus_path).exists():
                try:
                    import pandas as pd
                    from dash import dash_table

                    df = pd.read_csv(consensus_path, nrows=100)
                    return dash_table.DataTable(
                        data=df.to_dict("records"),
                        columns=[{"name": c, "id": c} for c in df.columns],
                        style_table={"overflowX": "auto"},
                        page_size=20,
                    )
                except Exception:
                    pass
            return html.Div("Taxonomy consensus table not available.", className="empty-state-sub")

        if stage == "phylo":
            phylo_out = outputs.get("phylo") or {}
            tree_path = (phylo_out.get("files") or {}).get("tree_nwk") or (phylo_out.get("files") or {}).get("newick")
            if tree_path and Path(tree_path).exists():
                rel = Path(tree_path).relative_to(runs_root) if _safe_relative(Path(tree_path), runs_root) else None
                return html.Div([
                    html.Div("Tree file available for download.", className="status-pill"),
                    html.A("Download .nwk", href=f"/runs/{rel.as_posix()}" if rel else "#"),
                ])
            return html.Div("Phylogenetic tree not available.", className="empty-state-sub")

        return html.Div(f"No preview available for {stage}.", className="empty-state-sub")

    @app.callback(
        Output("results-stats", "children"),
        Output("results-exports", "children"),
        Output("results-logs", "children"),
        Input("store-qc-outputs", "data"),
        Input("store-preprocess-outputs", "data"),
        Input("store-analysis-outputs", "data"),
        Input("store-databases-outputs", "data"),
        Input("store-phylo-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_results(qc_out, pre_out, analysis_out, db_out, phylo_out, run_config_data):
        qc = _normalize_store(qc_out)
        pre = _normalize_store(pre_out)
        analysis = _normalize_store(analysis_out)
        db = _normalize_store(db_out)
        phylo = _normalize_store(phylo_out)
        enabled = (run_config_data or {}).get("enabled_stages", {}) if isinstance(run_config_data, dict) else {}
        run_id = (run_config_data or {}).get("run_id", "")
        output_dir = (run_config_data or {}).get("output_dir", "")
        stage_records = {
            "qc": qc,
            "preprocess": pre,
            "analysis": analysis,
            "databases": db,
            "phylo": phylo,
        }
        active_records = [record for name, record in stage_records.items() if enabled.get(name, True)]

        # Stats
        status = "Complete" if all(r.get("status") == "success" for r in active_records) else "Incomplete"
        n_samples = len((run_config_data or {}).get("samples") or [])
        n_stages = sum(1 for r in active_records if r.get("status") and r.get("status") != "skipped")
        stats = [html.Div([html.Div(status, className="stat-value"), html.Div("Run status", className="stat-label")]), html.Div([html.Div(str(n_samples), className="stat-value"), html.Div("Samples", className="stat-label")]), html.Div([html.Div(str(n_stages), className="stat-value"), html.Div("Stages", className="stat-label")])]

        # Exports: collect files under runs dir and link to flask route
        runs_root = Path(get_runs_dir())
        files = []
        for out in active_records:
            files.extend(_collect_files_from_outputs(out.get("outputs", {})))

        links = []
        for f in files:
            try:
                p = Path(f)
                rel = p.relative_to(runs_root)
                href = f"/runs/{rel.as_posix()}"
                links.append(html.Div(html.A(p.name, href=href, target="_blank")))
            except ValueError:
                continue

        if not links:
            exports = html.Div([html.Div("No exported artifacts yet.", className="empty-state-sub")])
        else:
            run_config_raw = run_config_data or {}
            run_dir = run_config_raw.get("output_dir", "")
            if run_dir:
                zip_link = html.Div(
                    html.A("⬇ Download full run (.zip)", href=f"/download/run_zip?run_dir={run_dir}", className="btn-parafish"),
                    style={"marginTop": "12px"},
                )
                links.append(zip_link)
            exports = html.Div([html.Div("Downloads", className="card-title"), html.Div(links)])

        # Logs: aggregate simple textual logs
        logs = []
        for record in active_records:
            for line in (record.get("log_lines") or []):
                logs.append(f"[{line.get('ts')}] {line.get('level')}: {line.get('msg')}")

        seen = set()
        unique_logs = []
        for line in logs:
            parts = line.split("] ", 1)
            if len(parts) != 2 or ": " not in parts[1]:
                key = (line, None, None)
            else:
                ts_part = parts[0].lstrip("[")
                level_msg = parts[1].split(": ", 1)
                key = (ts_part, level_msg[0], level_msg[1] if len(level_msg) > 1 else "")
            if key not in seen:
                seen.add(key)
                unique_logs.append(line)

        log_block = html.Pre("\n".join(unique_logs)) if unique_logs else html.Div("No logs available.")

        return stats, exports, log_block
