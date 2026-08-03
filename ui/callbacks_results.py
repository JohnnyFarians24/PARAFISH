"""Populate the Results page from stored step outputs and expose download links."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, html
from dash.exceptions import PreventUpdate

from core.config import get_runs_dir
from ui.components.table import data_table
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


# Per-stage marker keys: a stage record's ``outputs`` mixes per-sample entries with
# run-wide artifacts (``vsearch_pool``, ``otu_matrix``, ``taxonomy``, plain paths...).
# A sample entry is a dict carrying that stage's per-sample step.
_SAMPLE_MARKERS = {
    "qc": ("metrics", "files"),
    "preprocess": ("trimmomatic", "pear"),
    "analysis": ("vsearch_filter",),
}

_NON_SAMPLE_KEYS = {"vsearch_pool", "dada2", "taxonomy", "phylo"}


def _sample_ids(stage: str, outputs: dict) -> list[str]:
    """Return the real per-sample keys of a stage record's outputs."""
    markers = _SAMPLE_MARKERS.get(stage)
    if not markers:
        return []
    return [
        key
        for key, value in (outputs or {}).items()
        if key not in _NON_SAMPLE_KEYS and isinstance(value, dict) and any(marker in value for marker in markers)
    ]


def _find_abundance_matrix(analysis_outputs: dict, run_output_dir) -> str | None:
    """Locate a feature x sample abundance matrix for cross-sample comparison.

    Prefers the VSEARCH ``otu_matrix`` (already CSV), then its raw TSV, then the
    DADA2 abundance table, then well-known paths under the run output directory.
    """
    for key in ("otu_matrix", "otu_table_tsv"):
        value = (analysis_outputs or {}).get(key)
        if isinstance(value, str) and value and Path(value).exists():
            return value
    dada2_files = ((analysis_outputs or {}).get("dada2") or {}).get("files") or {}
    abundance = dada2_files.get("asv_abundance_csv")
    if isinstance(abundance, str) and abundance and Path(abundance).exists():
        return abundance
    if run_output_dir:
        for rel in ("analysis/otu_matrix.csv", "analysis/pooled/otu_table.tsv", "analysis/dada2/asv_abundance.csv"):
            candidate = Path(run_output_dir) / rel
            if candidate.exists():
                return str(candidate)
    return None


def _run_link(path_value, label: str, runs_root: Path):
    """A download link for a file inside the runs root, or ``None``."""
    if not path_value:
        return None
    path = Path(str(path_value))
    if not path.exists():
        return None
    try:
        rel = path.relative_to(runs_root)
    except ValueError:
        return None
    return html.Div(html.A(label, href=f"/runs/{rel.as_posix()}", target="_blank"), className="download-row")


def _control_samples(samples_data) -> list[str]:
    """Internal ids of samples the user marked as negative controls."""
    raw = (samples_data or {}).get("samples") if isinstance(samples_data, dict) else samples_data
    return [str(sample.get("internal_id")) for sample in (raw or []) if isinstance(sample, dict) and sample.get("is_control")]


_GRAPH_CONFIG = {"displaylogo": False, "responsive": True}


def register_callbacks(app: Any) -> None:

    @app.callback(
        Output("results-sankey", "children"),
        Input("store-analysis-outputs", "data"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_sankey(analysis_out, db_out, run_config_data):
        from dash import dcc

        from analysis.results_summary import clustered_read_total, read_fate_sankey, species_classified_total
        from ui.plots import read_fate_sankey_figure

        output_dir = (run_config_data or {}).get("output_dir")
        if not output_dir:
            return html.Div("Run the pipeline to see the read-fate diagram.", className="empty-state-sub")

        # Both totals are read from the tables the stages already persisted. Deriving
        # them here would mean scanning clusters.uc (one row per dereplicated read —
        # millions on a real run) on every render, which is what made the server stop
        # responding.
        sankey = read_fate_sankey(
            output_dir,
            feature_reads=clustered_read_total(output_dir),
            classified_reads=species_classified_total(output_dir),
        )
        figure = read_fate_sankey_figure(sankey)
        if figure is None:
            return html.Div("Read-fate diagram needs the preprocess metrics; run the pipeline first.", className="empty-state-sub")
        return dcc.Graph(figure=figure, config=_GRAPH_CONFIG)

    @app.callback(
        Output("results-dashboard", "children"),
        Input("store-analysis-outputs", "data"),
        Input("store-databases-outputs", "data"),
        Input("store-samples", "data"),
        Input("results-target-input", "value"),
        State("store-run-config", "data"),
    )
    def _render_dashboard(analysis_out, db_out, samples_data, target_text, run_config_data):
        from dash import dash_table, dcc

        from analysis.results_summary import sample_dashboard
        from ui.plots import dashboard_richness_figure

        output_dir = (run_config_data or {}).get("output_dir")
        if not output_dir:
            return html.Div("Run the pipeline to see the sample dashboard.", className="empty-state-sub")
        targets = [taxon.strip() for taxon in (target_text or "").split(",") if taxon.strip()]
        rows = sample_dashboard(output_dir, control_samples=_control_samples(samples_data), target_taxa=targets)
        if not rows:
            return html.Div("No per-sample metrics yet.", className="empty-state-sub")

        children = []
        figure = dashboard_richness_figure(rows)
        if figure is not None:
            children.append(dcc.Graph(figure=figure, config=_GRAPH_CONFIG))
        children.append(
            data_table(
                data=[{**row, "control": "yes" if row["is_control"] else ""} for row in rows],
                columns=[
                    {"name": "Sample", "id": "sample"},
                    {"name": "Control", "id": "control"},
                    {"name": "Reads", "id": "reads", "type": "numeric", "format": {"specifier": ","}},
                    {"name": "Richness", "id": "richness", "type": "numeric"},
                    {"name": "Shannon", "id": "shannon", "type": "numeric"},
                    {"name": "% to species", "id": "pct_species", "type": "numeric"},
                    {"name": "Target taxa", "id": "n_targets", "type": "numeric"},
                ],
                sort_action="native",
                page_size=25,
                conditional=[{"if": {"filter_query": '{control} = "yes"'}, "backgroundColor": "rgba(248,81,73,0.10)"}],
            )
        )
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "10px"})

    @app.callback(
        Output("results-targets", "children"),
        Input("results-target-input", "value"),
        Input("store-databases-outputs", "data"),
        State("store-samples", "data"),
        State("store-run-config", "data"),
    )
    def _render_targets(target_text, db_out, samples_data, run_config_data):
        from dash import dash_table

        from analysis.results_summary import target_taxon_report

        targets = [taxon.strip() for taxon in (target_text or "").split(",") if taxon.strip()]
        if not targets:
            return html.Div("Enter one or more taxa above to build their report cards.", className="empty-state-sub")
        output_dir = (run_config_data or {}).get("output_dir")
        if not output_dir:
            return html.Div("Run the pipeline first.", className="empty-state-sub")
        rows = target_taxon_report(output_dir, targets, control_samples=_control_samples(samples_data))

        def _fmt_prev(row):
            if row["prevalence"] is None:
                return "not detected"
            ci = f" (95% CI {row['ci_low'] * 100:.0f}–{row['ci_high'] * 100:.0f}%)" if row["ci_low"] is not None else ""
            return f"{row['prevalence'] * 100:.0f}%{ci}"

        data = [
            {
                "taxon": row["taxon"],
                "detected_in": row["detected_in"],
                "prevalence": _fmt_prev(row),
                "n": f"{row['n_samples_detected']}/{row['n_samples']}" if row["n_samples"] else str(row["n_samples_detected"]),
                "total_reads": row["total_reads"],
                "best_reliability": row["best_reliability"],
            }
            for row in rows
        ]
        return data_table(
            data=data,
            columns=[
                {"name": "Taxon", "id": "taxon"},
                {"name": "Detected in", "id": "detected_in"},
                {"name": "Prevalence", "id": "prevalence"},
                {"name": "Samples", "id": "n"},
                {"name": "Total reads", "id": "total_reads", "type": "numeric", "format": {"specifier": ","}},
                {"name": "Best support", "id": "best_reliability"},
            ],
            conditional=[
                {"if": {"filter_query": '{best_reliability} = "high"', "column_id": "best_reliability"}, "color": "#3fb950", "fontWeight": "600"},
                {"if": {"filter_query": '{best_reliability} = "low"', "column_id": "best_reliability"}, "color": "#f85149", "fontWeight": "600"},
                {"if": {"filter_query": '{prevalence} = "not detected"'}, "color": "#8b98aa"},
            ],
        )

    @app.callback(
        Output("results-comparison", "children"),
        Input("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_comparison(analysis_out, run_config_data):
        from dash import dcc
        from analysis.sample_comparison import build_comparison
        from ui.plots import comparison_figures

        analysis = _normalize_store(analysis_out)
        outputs = analysis.get("outputs") or {}
        matrix_path = _find_abundance_matrix(outputs, (run_config_data or {}).get("output_dir"))
        if not matrix_path:
            return html.Div("No cross-sample abundance matrix available yet. Run the analysis stage to compare samples.", className="empty-state-sub")
        data = build_comparison(matrix_path)
        if data.get("error"):
            return html.Div(f"Comparison unavailable: {data['error']}", className="empty-state-sub")
        figures = comparison_figures(data)
        if not figures:
            return html.Div("Not enough data to compare samples.", className="empty-state-sub")
        header = html.Div(
            f"{data['n_samples']} sample(s) · {data['n_features']} features · {data['shared']} shared by all samples",
            className="mini-note",
        )
        graphs = [dcc.Graph(figure=fig, config={"displaylogo": False, "responsive": True}) for fig in figures]
        return html.Div([header, *graphs], style={"display": "flex", "flexDirection": "column", "gap": "12px"})

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
        Output("results-stage-selector", "value"),
        Output("results-stage-selector", "disabled"),
        Output("results-stage-selector", "placeholder"),
        Input("store-run-config", "data"),
        Input("store-project", "data"),
        State("results-stage-selector", "value"),
    )
    def _populate_stage_selector(run_config_data, project_data, current):
        """Stages come from the run on disk, and ALL of them are listed.

        Two bugs met here. The options were gated on ``status == "success"``, so a run
        whose analysis failed offered nothing to select and looked like a dead control
        rather than a failed run. And they were read from ``session`` stores that only
        the progress poller ever fills — so after a restart, in a new tab, or simply on
        arriving at Results later, the dropdown was empty no matter what had run.
        """
        from core.run_registry import STATUS_SUFFIX, load_run_state, run_dir_for, run_stages

        run_dir = run_dir_for(run_config_data, project_data)
        state = load_run_state(run_dir)
        if not state:
            return [], None, True, "No active run"

        stages = run_stages(state)
        options = [
            {
                # The status travels in the label, so a failed or skipped stage is
                # selectable AND self-explanatory.
                "label": f"{stage['label']}{STATUS_SUFFIX.get(stage['status'], '')}",
                "value": stage["key"],
                "disabled": not stage["has_outputs"],
            }
            for stage in stages
        ]
        selectable = [stage["key"] for stage in stages if stage["has_outputs"]]
        value = current if current in selectable else (selectable[0] if selectable else None)
        return options, value, False, "Choose a stage..."

    @app.callback(
        Output("results-sample-selector", "options"),
        Output("results-sample-selector", "value"),
        Output("results-sample-selector", "disabled"),
        Output("results-sample-selector", "placeholder"),
        Input("results-stage-selector", "value"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _populate_sample_selector(stage, run_config_data, project_data):
        """Chained off the stage: not every stage has per-sample output."""
        from core.run_registry import PER_SAMPLE_STAGES, list_run_samples, load_run_state, run_dir_for

        if not stage:
            return [], None, True, "Choose a stage first"

        state = load_run_state(run_dir_for(run_config_data, project_data))
        if not state:
            return [], None, True, "No active run"

        if stage not in PER_SAMPLE_STAGES:
            # Databases and Phylo run on the pooled feature set. Saying so beats an
            # enabled-but-empty dropdown the user keeps clicking at.
            return [], None, True, "This stage has no per-sample output"

        samples = list_run_samples(state, stage)
        if not samples:
            return [], None, True, "No samples in this stage"

        overview = {"label": "Overview (all samples)", "value": "__ALL__"}
        return [overview, *samples], "__ALL__", False, "Choose a sample..."

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

        if stage == "qc":
            # With no sample picked ("Overview") this used to fall through to "no
            # preview available", which is what the stage/sample selectors appeared
            # to break on. Overview lists every sample's report instead.
            if not sample_id:
                links = []
                for sid in _sample_ids("qc", outputs):
                    for path in ((outputs.get(sid) or {}).get("files") or {}).values():
                        if str(path).endswith(".html"):
                            link = _run_link(path, f"{sid}: {Path(str(path)).name}", runs_root)
                            if link is not None:
                                links.append(link)
                if not links:
                    return html.Div("No FastQC reports found.", className="empty-state-sub")
                return html.Div(
                    [html.Div("FastQC reports — pick a sample above to view one inline.", className="mini-note"), *links],
                    className="download-list",
                )
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

        if stage == "preprocess":
            if not sample_id:
                links = []
                for sid in _sample_ids("preprocess", outputs):
                    entry = outputs.get(sid) or {}
                    merged = ((entry.get("pear") or {}).get("files") or {}).get("assembled")
                    trimmed = ((entry.get("trimmomatic") or {}).get("files") or {}).get("r1_trimmed")
                    link = _run_link(merged or trimmed, f"{sid}: {Path(str(merged or trimmed)).name}" if (merged or trimmed) else sid, runs_root)
                    if link is not None:
                        links.append(link)
                if not links:
                    return html.Div("No preprocess outputs found.", className="empty-state-sub")
                return html.Div(
                    [html.Div("Merged (or trimmed) reads per sample. Pick a sample above for its full file list.", className="mini-note"), *links],
                    className="download-list",
                )
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
            if sample_id:
                # The per-sample step is stored under "vsearch_filter"; looking for a
                # "vsearch" key meant selecting a sample never showed anything.
                filter_out = (outputs.get(sample_id) or {}).get("vsearch_filter") or {}
                fasta_path = (filter_out.get("files") or {}).get("filtered_fasta")
                label = f"filtered reads for {sample_id}"
            else:
                fasta_path = outputs.get("combined_fasta")
                label = "OTUs/sequences"
            if fasta_path and Path(str(fasta_path)).exists():
                from core.utils import _count_fasta_seqs

                count = _count_fasta_seqs(Path(fasta_path))
                try:
                    rel = Path(str(fasta_path)).relative_to(runs_root)
                    href = f"/runs/{rel.as_posix()}"
                except ValueError:
                    href = "#"
                return html.Div([
                    html.Div(f"{count} {label}", className="stat-value"),
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
                    return data_table(
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

    # NOTE: results-exports is deliberately NOT an output here. The Results summary owns
    # it (ui/callbacks_results_summary.py), where it lists artifact sizes, a "download
    # all" zip, and — importantly — the artifacts that are MISSING, with the reason. Two
    # callbacks writing the same output is a hard Dash error ("Duplicate callback
    # outputs") and takes the whole page down.
    @app.callback(
        Output("results-stats", "children"),
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

        # Exports: collect files under runs dir, grouped by stage into collapsible
        # sections so the download hub reads as organised layers instead of one
        # long accumulated list.
        runs_root = Path(get_runs_dir())
        stage_labels = {
            "qc": "Quality control",
            "preprocess": "Preprocessing",
            "analysis": "Clustering / OTUs · ASVs",
            "databases": "Taxonomy & databases",
            "phylo": "Phylogenetics",
        }
        seen_paths: set[str] = set()
        stage_sections = []
        total_links = 0
        for stage_name, record in stage_records.items():
            if not enabled.get(stage_name, True):
                continue
            stage_links = []
            for f in _collect_files_from_outputs(record.get("outputs", {})):
                try:
                    p = Path(f)
                    key = str(p.resolve())
                    if key in seen_paths:
                        continue
                    seen_paths.add(key)
                    rel = p.relative_to(runs_root)
                    stage_links.append(html.Div(html.A(p.name, href=f"/runs/{rel.as_posix()}", target="_blank"), className="download-row"))
                except ValueError:
                    continue
            if stage_links:
                total_links += len(stage_links)
                stage_sections.append(
                    html.Details(
                        [
                            html.Summary(
                                [
                                    html.Span(stage_labels.get(stage_name, stage_name.title())),
                                    html.Span(f"{len(stage_links)}", className="download-count"),
                                ],
                                className="download-summary",
                            ),
                            html.Div(stage_links, className="download-list"),
                        ],
                        className="download-group",
                        open=(stage_name in ("databases", "phylo")),
                    )
                )

        if not total_links:
            exports = html.Div([html.Div("No exported artifacts yet.", className="empty-state-sub")])
        else:
            children = [html.Div("Downloads", className="card-title")]
            run_dir = (run_config_data or {}).get("output_dir", "")
            if run_dir:
                children.append(
                    html.Div(
                        html.A("⬇ Download full run (.zip)", href=f"/download/run_zip?run_dir={run_dir}", className="btn-parafish"),
                        style={"margin": "6px 0 14px"},
                    )
                )
            children.extend(stage_sections)
            exports = html.Div(children)

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

        del exports  # rendered by the Results summary, which owns that output
        return stats, log_block
