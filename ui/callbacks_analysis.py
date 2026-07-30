"""Callbacks related to the analysis stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, dash_table, dcc, html, no_update
from dash.exceptions import PreventUpdate

from analysis.analysis_metrics import DEFAULT_CORE_FRACTION, build_metrics_from_outputs
from analysis.dada2_runner import run_dada2
from analysis.vsearch_runner import run_vsearch_filter_sample, run_vsearch_pool_and_cluster
from ui.components.table import data_table
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import run_config_from_store, samples_from_store, serialize_step_output, build_sample_from_record, latest_run_dir



_GRAPH_CONFIG = {"displaylogo": False, "responsive": True}

_EMPTY = "Run the analysis stage to see OTU/ASV outputs."


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





def _resolve_output_dir(run_config_data, project_data=None) -> str | None:
    """The run's output directory, from the session store or reconstructed from disk.

    ``run_config_data`` is per-session: a fresh tab or a restarted app has it empty, so
    ``output_dir`` is gone even though the run is on disk. Fall back to the project's newest
    run (``latest_run_dir``, ordered by run_id, not mtime) instead of blanking. Eixo B only:
    it resolves WHICH run, never which feature set — the ``otu`` assumption below is a
    separate (Eixo A) concern, tracked in the living registry, unblocked by an ASV golden.
    """
    output_dir = (run_config_data or {}).get("output_dir")
    if output_dir:
        return str(output_dir)
    run_dir = latest_run_dir(project_data) if project_data else None
    return str(run_dir) if run_dir else None


def _analysis_metric_paths(analysis_data, run_config_data, project_data=None) -> dict[str, str]:
    """Locate the analysis tables the runner persisted.

    Unlike QC and Preprocess there is no cheap rebuild here: rarefaction is a
    Monte-Carlo subsampling that belongs in the stage, not in a render. Runs
    predating this layer simply show the empty state.
    """
    record = _normalize_store(analysis_data)
    recorded = record.get("metrics_files") or {}
    paths = {key: value for key, value in recorded.items() if isinstance(value, str) and Path(value).exists()}
    if paths:
        return paths
    output_dir = _resolve_output_dir(run_config_data, project_data)
    if not output_dir:
        return {}
    from analysis.analysis_metrics import analysis_metrics_paths
    from ui.feature_sets import resolve_feature_set

    # These callbacks carry no set selector, so render whichever set the run produced —
    # asv for a DADA2-only run, otu otherwise. resolve_feature_set reads the recorded sets
    # (or falls to disk when the store is empty), so this stays right in a fresh session.
    feature_set = resolve_feature_set(None, record.get("outputs") or {}, output_dir)
    return {key: str(path) for key, path in analysis_metrics_paths(Path(output_dir) / "analysis", feature_set).items()}


def _download_link(path_value, label: str):
    """Render a /runs/... download link, or None when the file is outside the runs root."""
    from core.config import get_runs_dir

    if not path_value:
        return None
    path = Path(str(path_value))
    if not path.exists():
        return None
    try:
        rel = path.relative_to(Path(get_runs_dir()))
    except ValueError:
        return None
    return html.Div(html.A(f"{label}: {path.name}", href=f"/runs/{rel.as_posix()}", target="_blank"), className="download-row")


def _analysis_stat_values(analysis_data, run_config_data, databases_data=None, project_data=None) -> dict:
    """The three Analysis KPI values (OTUs / ASVs / Consensus) as raw strings + hints.

    Split out of the callback so a restored session can be tested directly — the OTU/ASV
    counts come from the stage record in ``store-analysis-outputs``, which ``restore_run_state``
    repopulates from disk with the run's file paths. A count that reads an absolute FASTA
    path the store carries is store-only but survives restore, because the store is rebuilt
    faithfully; the test locks that chain so a future change that stops repopulating the
    record fails loudly instead of silently blanking the cards.

    The Consensus card reads the DATABASES store, not this one: the consensus is a
    Databases-stage output and never reaches ``store-analysis-outputs``, so the old read of
    ``analysis_outputs.taxonomy`` could never fire its own ``"ready"`` branch — a crossed
    wire that made the card looks-alive-does-nothing. Resolved through the canonical
    ``_find_consensus_csv`` (per-feature-set → disk → flat pre-feature-set shape), across
    BOTH feature sets, so it never lights for OTU while staying dead for ASV.
    """
    from core.utils import _count_fasta_seqs

    data = _normalize_store(analysis_data)
    outputs = data.get("outputs") or {}
    method = str((run_config_data or {}).get("analysis_method") or "vsearch").lower()

    combined_fasta = outputs.get("combined_fasta") or ""
    vsearch_fasta = outputs.get("combined_fasta_vsearch") or ""
    otu_count = "—"
    for candidate in (vsearch_fasta, combined_fasta):
        if candidate and Path(candidate).exists():
            otu_count = _count_fasta_seqs(Path(candidate))
            break
    otu_hint = "" if method in ("vsearch", "both") else "VSEARCH not selected"

    dada2_out = outputs.get("dada2") or {}
    asv_fasta = (dada2_out.get("files") or {}).get("asv_sequences_fasta") or ""
    if asv_fasta and Path(asv_fasta).exists():
        asv_value, asv_hint = _count_fasta_seqs(Path(asv_fasta)), ""
    elif method in ("dada2", "both"):
        asv_value, asv_hint = "—", "DADA2 pending / failed"
    else:
        asv_value, asv_hint = "off", "DADA2 not selected"

    from ui.callbacks_databases import _find_consensus_csv

    db_outputs = _normalize_store(databases_data).get("outputs") or {}
    run_output_dir = _resolve_output_dir(run_config_data, project_data)
    if any(_find_consensus_csv(db_outputs, run_output_dir, feature_set) for feature_set in ("otu", "asv")):
        consensus_value, consensus_hint = "ready", ""
    else:
        consensus_value, consensus_hint = "—", "runs in Databases stage"

    return {
        "otu": otu_count, "otu_hint": otu_hint,
        "asv": asv_value, "asv_hint": asv_hint,
        "consensus": consensus_value, "consensus_hint": consensus_hint,
    }


def register_callbacks(app: Any) -> None:

    @app.callback(
        Output("analysis-occupancy", "children"),
        Input("store-analysis-outputs", "data"),
        Input("analysis-occupancy-set", "value"),
        Input("analysis-occupancy-filter", "value"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_occupancy(analysis_data, feature_set, mode, run_config_data, project_data):
        """Which samples each feature was found in — the cross-sample question."""
        import pandas as pd

        from analysis.analysis_metrics import analysis_metrics_paths
        from ui.feature_sets import FEATURE_SET_NOUN

        feature_set = feature_set or "otu"
        noun = FEATURE_SET_NOUN.get(feature_set, "feature")
        output_dir = _resolve_output_dir(run_config_data, project_data)
        if not output_dir:
            return html.Div(_EMPTY, className="empty-state-sub")

        paths = analysis_metrics_paths(Path(output_dir) / "analysis", feature_set)
        occupancy = paths.get("occupancy_csv")
        if not occupancy or not Path(occupancy).exists():
            return html.Div(
                f"Sem tabela de ocupação para {noun}s. Corre a análise com o método que gera {noun}s.",
                className="empty-state-sub",
            )
        try:
            frame = pd.read_csv(occupancy)
        except Exception as exc:
            return html.Div(f"Tabela indisponível: {exc}", className="empty-state-sub")
        if frame.empty:
            return html.Div("Sem features.", className="empty-state-sub")

        total = len(frame)
        shared = int(frame["is_shared"].sum()) if "is_shared" in frame.columns else 0
        if mode == "shared" and "is_shared" in frame.columns:
            frame = frame[frame["is_shared"] == 1]
        elif mode == "unique" and "is_unique" in frame.columns:
            frame = frame[frame["is_unique"] == 1]

        headers = {
            "feature": noun,
            "n_samples_present": "Nº amostras",
            "samples_present": "Amostras",
            "total_reads": "Reads",
            "max_sample": "Amostra dominante",
            "max_sample_reads": "Reads aí",
        }
        display = frame[[c for c in headers if c in frame.columns]]
        if display.empty:
            return html.Div("Nenhuma feature corresponde a este filtro.", className="empty-state-sub")

        return html.Div(
            [
                html.Div(
                    [
                        html.Div(f"{total} {noun}s no total", className="status-pill"),
                        html.Div(f"{shared} em mais de uma amostra", className="status-pill"),
                        html.Div(f"{total - shared} exclusivas de uma amostra", className="status-pill"),
                    ],
                    className="chip-list",
                ),
                data_table(
                    data=display.to_dict("records"),
                    columns=[{"name": headers.get(str(c), str(c)), "id": str(c)} for c in display.columns],
                    sort_action="native",
                    filter_action="native",
                    page_size=15,
                ),
            ],
            style={"display": "flex", "flexDirection": "column", "gap": "10px"},
        )

    @app.callback(
        Output("analysis-summary", "children", allow_duplicate=True),
        Input("store-analysis-outputs", "data"),
        State("store-samples", "data"),
        State("store-preprocess-outputs", "data"),
        prevent_initial_call="initial_duplicate",
    )
    def _sync_analysis_summary(analysis_data, samples_data, preprocess_data):
        # A pipeline launched from Setup fills the store directly, so the summary
        # has to follow the store and not only the on-page "Run Analysis" button.
        return run_analysis(
            samples_data=samples_data,
            step_outputs={"preprocess": _normalize_store(preprocess_data), "analysis": _normalize_store(analysis_data)},
        )

    @app.callback(
        Output("analysis-chimera", "children"),
        Input("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_chimera_stats(analysis_data, run_config_data, project_data):
        import pandas as pd

        paths = _analysis_metric_paths(analysis_data, run_config_data, project_data)
        if "chimera_stats_csv" not in paths:
            return html.Div("Run the analysis stage to see chimera detection results.", className="empty-state-sub")
        try:
            stats = pd.read_csv(paths["chimera_stats_csv"]).iloc[0]
        except Exception as exc:
            return html.Div(f"Chimera stats unavailable: {exc}", className="empty-state-sub")

        def _pill(value, label):
            return html.Div([html.Div(value, className="stat-value"), html.Div(label, className="stat-label")], className="stat-card")

        basis = str(stats.get("read_basis", "centroid"))
        note = (
            "Read shares use the cluster map, so the denominator is every clustered read."
            if basis == "cluster"
            else "No cluster map found: read shares fall back to the centroids' own dereplicated abundance."
        )
        return html.Div(
            [
                html.Div(
                    [
                        _pill(f"{stats['pct_chimeric_variants']:.2f}%", f"of variants ({int(stats['chimeric_variants'])}/{int(stats['total_variants'])})"),
                        _pill(f"{stats['pct_chimeric_reads']:.2f}%", f"of reads ({int(stats['chimeric_reads']):,}/{int(stats['total_reads']):,})"),
                    ],
                    className="grid-2",
                ),
                html.Div(note, className="mini-note"),
            ]
        )

    @app.callback(
        Output("analysis-dada2-tracking", "children"),
        Input("store-analysis-outputs", "data"),
    )
    def _render_dada2_tracking(analysis_data):
        import pandas as pd

        outputs = _normalize_store(analysis_data).get("outputs") or {}
        tracking = ((outputs.get("dada2") or {}).get("files") or {}).get("read_tracking_csv")
        if not tracking or not Path(str(tracking)).exists():
            return html.Div("DADA2 was not run, so there is no denoising read tracking.", className="empty-state-sub")
        try:
            frame = pd.read_csv(tracking)
        except Exception as exc:
            return html.Div(f"DADA2 read tracking unavailable: {exc}", className="empty-state-sub")
        if frame.empty:
            return html.Div("DADA2 produced no read tracking.", className="empty-state-sub")
        return html.Div(
            [
                html.Div("Reads surviving each DADA2 step, per sample.", className="mini-note"),
                data_table(
                    data=frame.to_dict("records"),
                    columns=[{"name": str(column), "id": str(column)} for column in frame.columns],
                    sort_action="native",
                    page_size=15,
                ),
            ]
        )

    @app.callback(
        Output("analysis-rarefaction-note", "children"),
        Output("analysis-alpha-table", "children"),
        Input("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_alpha_diversity(analysis_data, run_config_data, project_data):
        import pandas as pd

        paths = _analysis_metric_paths(analysis_data, run_config_data, project_data)
        if "alpha_diversity_csv" not in paths:
            return "", html.Div(_EMPTY, className="empty-state-sub")
        try:
            alpha = pd.read_csv(paths["alpha_diversity_csv"])
            depth_row = pd.read_csv(paths["rarefaction_depth_csv"]).iloc[0] if "rarefaction_depth_csv" in paths else None
        except Exception as exc:
            return "", html.Div(f"Alpha diversity unavailable: {exc}", className="empty-state-sub")

        note = ""
        if depth_row is not None:
            dropped = int(depth_row["n_samples_dropped"])
            note = (
                f"All samples rarefied to {int(depth_row['depth']):,} reads — richness grows with sequencing "
                f"effort, so alpha diversity is only comparable at equal depth."
            )
            if dropped:
                note += f" {dropped} empty sample(s) excluded."
        columns = [
            {"name": "Sample", "id": "sample"},
            {"name": "Reads", "id": "reads", "type": "numeric", "format": {"specifier": ","}},
            {"name": "Observed", "id": "observed", "type": "numeric"},
            {"name": "Chao1", "id": "chao1", "type": "numeric"},
            {"name": "Shannon", "id": "shannon", "type": "numeric"},
            {"name": "Simpson", "id": "simpson", "type": "numeric"},
            {"name": "Pielou", "id": "pielou", "type": "numeric"},
            {"name": "Hill q=0", "id": "hill_q0", "type": "numeric"},
            {"name": "Hill q=1", "id": "hill_q1", "type": "numeric"},
            {"name": "Hill q=2", "id": "hill_q2", "type": "numeric"},
        ]
        table = data_table(
            data=alpha.to_dict("records"),
            columns=columns,
            sort_action="native",
            page_size=15,
        )
        return note, table

    @app.callback(
        Output("analysis-upset", "children"),
        Input("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_upset(analysis_data, run_config_data, project_data):
        import pandas as pd

        from ui.plots import upset_figure

        paths = _analysis_metric_paths(analysis_data, run_config_data, project_data)
        if "intersections_csv" not in paths:
            return html.Div(_EMPTY, className="empty-state-sub")
        try:
            intersections = pd.read_csv(paths["intersections_csv"])
        except Exception as exc:
            return html.Div(f"Intersections unavailable: {exc}", className="empty-state-sub")
        samples = sorted({sample for combination in intersections["combination"] for sample in str(combination).split("|")})
        figure = upset_figure(intersections, samples)
        if figure is None:
            return html.Div("Not enough samples to compare feature sets.", className="empty-state-sub")

        children = [dcc.Graph(figure=figure, config=_GRAPH_CONFIG)]
        if "prevalence_csv" in paths:
            try:
                prevalence = pd.read_csv(paths["prevalence_csv"])
                n_core = int(prevalence["is_core"].sum())
                children.append(
                    html.Div(
                        f"{n_core} core feature(s) — present in at least "
                        f"{int(DEFAULT_CORE_FRACTION * 100)}% of samples, out of {len(prevalence)} total.",
                        className="mini-note",
                    )
                )
            except Exception:
                pass
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "8px"})

    @app.callback(
        Output("analysis-sample-table", "children"),
        Input("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_analysis_sample_table(analysis_data, run_config_data, project_data):
        from analysis.sample_comparison import build_comparison
        from ui.callbacks_results import _find_abundance_matrix

        outputs = _normalize_store(analysis_data).get("outputs") or {}
        # Per-sample reads/richness come from the abundance matrix, the only place
        # the runners actually record them (the filter step emits no metrics).
        matrix_path = _find_abundance_matrix(outputs, _resolve_output_dir(run_config_data, project_data))
        if not matrix_path:
            return html.Div(_EMPTY, className="empty-state-sub")
        data = build_comparison(matrix_path)
        if data.get("error"):
            return html.Div(f"Per-sample summary unavailable: {data['error']}", className="empty-state-sub")

        exclusive = data.get("exclusive") or {}
        rows = []
        for entry in data.get("diversity") or []:
            sample_id = entry["sample"]
            filter_step = (outputs.get(sample_id) or {}).get("vsearch_filter") or {}
            rows.append(
                {
                    "Sample": sample_id,
                    "Reads": entry["reads"],
                    "Richness (OTUs/ASVs)": entry["richness"],
                    "Shannon": entry["shannon"],
                    "Exclusive features": exclusive.get(sample_id, 0),
                    "Filter": "success" if filter_step.get("success", True) else "failure",
                }
            )
        if not rows:
            return html.Div(_EMPTY, className="empty-state-sub")
        return data_table(
            data=sorted(rows, key=lambda row: row["Reads"], reverse=True),
            columns=[{"name": c, "id": c} for c in ("Sample", "Reads", "Richness (OTUs/ASVs)", "Shannon", "Exclusive features", "Filter")],
            sort_action="native",
            page_size=15,
        )

    @app.callback(
        Output("analysis-otu-table", "children"),
        Input("store-analysis-outputs", "data"),
        Input("analysis-abundance-mode", "value"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_otu_table(analysis_data, mode, run_config_data, project_data):
        import pandas as pd

        from ui.callbacks_results import _find_abundance_matrix

        outputs = _normalize_store(analysis_data).get("outputs") or {}
        matrix_path = _find_abundance_matrix(outputs, _resolve_output_dir(run_config_data, project_data))
        if not matrix_path:
            return html.Div(_EMPTY, className="empty-state-sub")
        try:
            sep = "\t" if str(matrix_path).endswith((".tsv", ".txt")) else ","
            frame = pd.read_csv(matrix_path, sep=sep)
        except Exception as exc:
            return html.Div(f"Abundance table unavailable: {exc}", className="empty-state-sub")

        feature_column = frame.columns[0]
        sample_columns = list(frame.columns[1:])
        if mode == "relative" and sample_columns:
            # Column-wise: each sample's counts become its own composition, which is
            # what makes samples of different depth comparable.
            totals = frame[sample_columns].sum(axis=0).replace(0, pd.NA)
            frame = frame.assign(**{column: (100.0 * frame[column] / totals[column]).round(4) for column in sample_columns})

        max_rows = 500
        truncated = len(frame) > max_rows
        unit = "relative abundance (% of each sample's reads)" if mode == "relative" else "absolute read counts"
        note = f"{len(frame)} features x {len(sample_columns)} samples, {unit}."
        if truncated:
            note += f" Showing the first {max_rows}."
        return html.Div(
            [
                html.Div(note, className="mini-note"),
                data_table(
                    data=frame.head(max_rows).to_dict("records"),
                    columns=[{"name": str(feature_column), "id": str(feature_column)}]
                    + [{"name": str(c), "id": str(c), "type": "numeric"} for c in sample_columns],
                    sort_action="native",
                    filter_action="native",
                    page_size=15,
                ),
            ]
        )

    @app.callback(
        Output("analysis-plots", "children"),
        Input("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_analysis_plots(analysis_data, run_config_data, project_data):
        import pandas as pd

        from analysis.sample_comparison import build_comparison
        from ui.callbacks_results import _find_abundance_matrix
        from ui.plots import (
            alpha_diversity_figure,
            comparison_figures,
            feature_length_figure,
            rarefaction_figure,
            rank_abundance_figure,
        )

        paths = _analysis_metric_paths(analysis_data, run_config_data, project_data)
        figures = []
        try:
            if "rarefaction_csv" in paths:
                depth = None
                if "rarefaction_depth_csv" in paths:
                    depth = int(pd.read_csv(paths["rarefaction_depth_csv"]).iloc[0]["depth"])
                figures.append(rarefaction_figure(pd.read_csv(paths["rarefaction_csv"]), depth))
            if "alpha_diversity_csv" in paths:
                figures.append(alpha_diversity_figure(pd.read_csv(paths["alpha_diversity_csv"])))
            if "rank_abundance_csv" in paths:
                figures.append(rank_abundance_figure(pd.read_csv(paths["rank_abundance_csv"])))
            if "feature_length_csv" in paths:
                figures.append(feature_length_figure(pd.read_csv(paths["feature_length_csv"])))
        except Exception as exc:
            return html.Div(f"Plots unavailable: {exc}", className="empty-state-sub")

        header = None
        outputs = _normalize_store(analysis_data).get("outputs") or {}
        matrix_path = _find_abundance_matrix(outputs, _resolve_output_dir(run_config_data, project_data))
        if matrix_path:
            data = build_comparison(matrix_path)
            if not data.get("error"):
                header = html.Div(
                    f"{data['n_samples']} sample(s) · {data['n_features']} features · {data['shared']} shared by all samples",
                    className="mini-note",
                )
                figures.extend(comparison_figures(data))

        figures = [figure for figure in figures if figure is not None]
        if not figures:
            return html.Div(_EMPTY, className="empty-state-sub")
        graphs = [dcc.Graph(figure=fig, config=_GRAPH_CONFIG) for fig in figures]
        children = [header, *graphs] if header is not None else graphs
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "12px"})

    @app.callback(
        Output("analysis-downloads", "children"),
        Input("store-analysis-outputs", "data"),
    )
    def _render_analysis_downloads(analysis_data):
        outputs = _normalize_store(analysis_data).get("outputs") or {}
        dada2_files = (outputs.get("dada2") or {}).get("files") or {}
        pool_files = (outputs.get("vsearch_pool") or {}).get("files") or {}
        wanted = [
            (outputs.get("combined_fasta_vsearch"), "OTU representatives (FASTA)"),
            (outputs.get("otu_table_tsv"), "OTU table (TSV)"),
            (outputs.get("otu_matrix"), "OTU matrix (CSV)"),
            (pool_files.get("otu_id_map_csv"), "OTU id → representative read (CSV)"),
            (outputs.get("combined_otu_uc"), "Cluster map (UC)"),
            (dada2_files.get("asv_sequences_fasta"), "ASV sequences (FASTA)"),
            (dada2_files.get("asv_abundance_csv"), "ASV abundance (CSV)"),
        ]
        links = [link for link in (_download_link(path, label) for path, label in wanted) if link is not None]
        if not links:
            return html.Div("No analysis artifacts yet.", className="empty-state-sub")
        return html.Div(links, className="download-list")

    @app.callback(
        Output("analysis-stat-otus", "children"),
        Output("analysis-stat-asvs", "children"),
        Output("analysis-stat-consensus", "children"),
        Input("store-analysis-outputs", "data"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _update_analysis_stats(analysis_data, databases_data, run_config_data, project_data):
        def _card(value, label, hint=""):
            children = [html.Div(str(value), className="stat-value"), html.Div(label, className="stat-label")]
            if hint:
                children.append(html.Div(hint, className="stat-hint"))
            return children

        values = _analysis_stat_values(analysis_data, run_config_data, databases_data, project_data)
        return (
            _card(values["otu"], "OTUs", values["otu_hint"]),
            _card(values["asv"], "ASVs", values["asv_hint"]),
            _card(values["consensus"], "Consensus", values["consensus_hint"]),
        )

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
        for sample_record in (preprocess_record.get("next_samples") or []):
            sample = build_sample_from_record(sample_record)
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

        run_vsearch_branch = run_config.analysis_method in ("vsearch", "both")
        run_dada2_branch = run_config.analysis_method in ("dada2", "both")

        log_lines = [{"ts": utc_now_iso(), "level": "RUN", "msg": f"Analysis method: {run_config.analysis_method}"}]
        outputs = {}
        success = True

        if run_vsearch_branch:
            filtered_fastas = {}
            for sample_record in (preprocess_record.get("next_samples") or []):
                sample_id = sample_record["internal_id"]
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
                filter_res = run_vsearch_filter_sample(assembled_path, run_config, run_config.output_dir / "analysis" / sample_id, sample_id)
                outputs.setdefault(sample_id, {})["vsearch_filter"] = serialize_step_output(filter_res)
                if filter_res.success:
                    filtered_fastas[sample_id] = filter_res.files["filtered_fasta"]
                    log_lines.append({"ts": utc_now_iso(), "level": "OK", "msg": f"VSEARCH filter completed: {sample_id}"})
                else:
                    success = False
                    log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"VSEARCH filter failed: {sample_id} -> {filter_res.error_message or 'unknown error'}"})

            if filtered_fastas:
                pool_res = run_vsearch_pool_and_cluster(filtered_fastas, run_config, run_config.output_dir / "analysis" / "pooled")
                outputs["vsearch_pool"] = serialize_step_output(pool_res)
                if pool_res.success:
                    outputs["combined_fasta_vsearch"] = str(pool_res.files["otus_fasta"])
                    log_lines.append({"ts": utc_now_iso(), "level": "OK", "msg": "VSEARCH pooled clustering completed"})
                else:
                    success = False
                    outputs["combined_fasta_vsearch"] = ""
                    log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"VSEARCH pooled clustering failed -> {pool_res.error_message or 'unknown error'}"})
            else:
                success = False
                outputs["combined_fasta_vsearch"] = ""
                log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": "Nenhuma amostra produziu reads filtradas para VSEARCH."})

        if run_dada2_branch:
            analysis_samples = [build_sample_from_record(sample_record) for sample_record in (preprocess_record.get("next_samples") or [])]
            dada2_result = run_dada2(analysis_samples, run_config, run_config.output_dir / "analysis" / "dada2")
            outputs["dada2"] = serialize_step_output(dada2_result)
            if dada2_result.success:
                outputs["combined_fasta_dada2"] = str(dada2_result.files.get("asv_sequences_fasta"))
                log_lines.append({"ts": utc_now_iso(), "level": "OK", "msg": "DADA2 completed"})
            else:
                success = False
                outputs["combined_fasta_dada2"] = ""
                log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"DADA2 failed -> {dada2_result.error_message or 'unknown error'}"})

        if run_vsearch_branch:
            outputs["combined_fasta"] = outputs.get("combined_fasta_vsearch") or ""
        elif run_dada2_branch:
            outputs["combined_fasta"] = outputs.get("combined_fasta_dada2") or ""
        else:
            outputs["combined_fasta"] = ""
            success = False
            log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"Analysis method inválido: {run_config.analysis_method}"})

        try:
            metrics_files = build_metrics_from_outputs(run_config.output_dir / "analysis", outputs)
        except Exception as exc:
            metrics_files = {}
            log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"Falha ao construir métricas de análise: {exc}"})

        record = {
            "signature": signature,
            "status": "success" if success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": outputs,
            "metrics_files": metrics_files,
            "log_lines": log_lines,
        }
        return run_analysis(samples_data=samples_data, step_outputs={"preprocess": preprocess_record, "analysis": record}), record
