"""Callbacks related to the phylogenetics stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio import SeqIO
from dash import Input, Output, State, dcc, html, no_update
from dash.exceptions import PreventUpdate

from phylo import phylo_input_error
from phylo.phylo_metrics import build_metrics_from_files
from phylo.phylo_runner import run_phylo as run_phylo_step
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import run_config_from_store


_GRAPH_CONFIG = {"displaylogo": False, "responsive": True}


def _phylo_files(phylo_data) -> dict:
    """The phylo stage's serialized file map."""
    data = _normalize_store(phylo_data)
    return ((data.get("outputs") or {}).get("phylo") or {}).get("files") or {}


def _phylo_metric_paths(phylo_data, run_config_data) -> dict[str, str]:
    """Locate the phylo tables the runner persisted."""
    record = _normalize_store(phylo_data)
    recorded = record.get("metrics_files") or {}
    paths = {key: value for key, value in recorded.items() if isinstance(value, str) and Path(value).exists()}
    if paths:
        return paths
    output_dir = (run_config_data or {}).get("output_dir")
    if not output_dir:
        return {}
    from phylo.phylo_metrics import phylo_metrics_paths
    from ui.feature_sets import phylo_dir

    # The metrics live under phylo/<set>/; fall back to the flat layout for older runs.
    for feature_set in ("otu", "asv"):
        base = phylo_dir(output_dir, feature_set)
        resolved = {key: str(path) for key, path in phylo_metrics_paths(base).items()}
        if resolved:
            return resolved
    return {}


def _mapping_csv(files: dict, run_config_data) -> str | None:
    """Path to ``reversed_asvs.csv``, which maps tree tips back to feature ids."""
    value = files.get("reversed_asvs_csv")
    if value and Path(str(value)).exists():
        return str(value)
    output_dir = (run_config_data or {}).get("output_dir")
    if output_dir:
        from ui.feature_sets import phylo_artifact

        for feature_set in ("otu", "asv"):
            candidate = phylo_artifact(output_dir, feature_set, "reversed_asvs_csv")
            if candidate:
                return candidate
    return None


def _tip_abundance(mapping: dict, analysis_data, run_config_data) -> dict:
    """Total reads behind each tree tip, keyed by tip name."""
    import pandas as pd

    from ui.callbacks_results import _find_abundance_matrix

    if not mapping:
        return {}
    analysis = _normalize_store(analysis_data)
    matrix_path = _find_abundance_matrix(analysis.get("outputs") or {}, (run_config_data or {}).get("output_dir"))
    if not matrix_path:
        return {}
    try:
        separator = "\t" if str(matrix_path).endswith((".tsv", ".txt")) else ","
        frame = pd.read_csv(matrix_path, sep=separator)
    except Exception:
        return {}
    feature_column, *sample_columns = frame.columns
    totals: dict[str, int] = {}
    for _, row in frame.iterrows():
        tip = mapping.get(str(row[feature_column]).split(";")[0].strip())
        if tip:
            totals[tip] = totals.get(tip, 0) + int(row[sample_columns].sum())
    return totals


def _tip_taxonomy(mapping: dict, databases_data, run_config_data) -> dict:
    """Best available taxon label for each tree tip, keyed by tip name."""
    import pandas as pd

    from ui.callbacks_databases import _find_consensus_csv

    if not mapping:
        return {}
    databases = _normalize_store(databases_data)
    consensus_path = _find_consensus_csv(databases.get("outputs") or {}, (run_config_data or {}).get("output_dir"))
    if not consensus_path:
        return {}
    try:
        consensus = pd.read_csv(consensus_path)
    except Exception:
        return {}
    ranks = ["Consensus_Species", "Consensus_Genus", "Consensus_Family", "Consensus_Order", "Consensus_Phylum"]
    labels: dict[str, str] = {}
    for _, row in consensus.iterrows():
        tip = mapping.get(str(row.get("ASV_ID", "")).split(";")[0].strip())
        if not tip:
            continue
        for rank in ranks:
            value = str(row.get(rank, "") or "").strip()
            if value and value.lower() not in {"nan", "unclassified"}:
                labels[tip] = value
                break
    return labels


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

        guard_error = phylo_input_error(seq_count)
        if guard_error:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": guard_error}],
                "error_message": guard_error,
            }
            return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

        # Write per feature set, exactly as the pipeline does — a flat phylo/ dir made
        # the per-set readers (tree viewer, KPIs, Results) look one directory too high.
        from ui.feature_sets import resolve_feature_set

        feature_set = resolve_feature_set(None, analysis_record.get("outputs") or {}, run_config.output_dir)
        phylo_base = run_config.output_dir / "phylo" / feature_set
        phylo_result = run_phylo_step(fasta_path, run_config, phylo_base)
        phylo_files = {key: str(value) for key, value in phylo_result.files.items()}
        log_lines = [{"ts": utc_now_iso(), "level": "OK" if phylo_result.success else "ERR", "msg": phylo_result.error_message or "Phylo pipeline completed"}]
        try:
            from ui.callbacks_results import _find_abundance_matrix

            metrics_files = build_metrics_from_files(
                phylo_base,
                phylo_files,
                _find_abundance_matrix(analysis_record.get("outputs") or {}, run_config.output_dir),
            )
        except Exception as exc:
            metrics_files = {}
            log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"Falha ao construir métricas de filogenia: {exc}"})

        record = {
            "signature": signature,
            "status": "success" if phylo_result.success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": {"phylo": {"files": phylo_files, "metrics": phylo_result.metrics}},
            "metrics_files": metrics_files,
            "log_lines": log_lines,
            "error_message": phylo_result.error_message,
        }
        return _build_phylo_summary(samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": databases_record, "phylo": record}), record

    @app.callback(
        Output("phylo-alignment", "children"),
        Input("store-phylo-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_alignment_qc(phylo_data, run_config_data):
        import pandas as pd

        from ui.plots import alignment_profile_figure

        paths = _phylo_metric_paths(phylo_data, run_config_data)
        if "alignment_summary_csv" not in paths:
            return html.Div("Run the phylo stage to see alignment QC.", className="empty-state-sub")
        try:
            summary = pd.read_csv(paths["alignment_summary_csv"]).iloc[0]
        except Exception as exc:
            return html.Div(f"Alignment QC unavailable: {exc}", className="empty-state-sub")

        def _pill(value, label, color=None):
            return html.Div(
                [html.Div(str(value), className="stat-value", style={"color": color} if color else {}), html.Div(label, className="stat-label")],
                className="stat-card",
            )

        n_reversed = int(summary.get("n_reversed", 0) or 0)
        children = [
            html.Div(
                [
                    _pill(int(summary["n_sequences"]), "Sequences"),
                    _pill(f"{int(summary['alignment_length']):,}", "Alignment columns"),
                    _pill(f"{summary['pct_gaps']:.1f}%", "Gaps", color="var(--amber)" if summary["pct_gaps"] > 50 else None),
                    _pill(n_reversed, "Reverse-complemented", color="var(--amber)" if n_reversed else None),
                ],
                className="grid-4",
            )
        ]
        if n_reversed:
            children.append(
                html.Div(
                    f"MAFFT reverse-complemented {n_reversed} sequence(s): {summary.get('reversed_ids', '')}. "
                    "Mixed orientation in the input is worth checking before trusting the tree.",
                    className="mini-note",
                    style={"color": "var(--amber)"},
                )
            )
        if "alignment_columns_csv" in paths:
            try:
                figure = alignment_profile_figure(pd.read_csv(paths["alignment_columns_csv"]))
            except Exception:
                figure = None
            if figure is not None:
                children.append(dcc.Graph(figure=figure, config=_GRAPH_CONFIG))
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "10px"})

    @app.callback(
        Output("phylo-tree-viewer", "children"),
        Input("store-phylo-outputs", "data"),
        Input("store-analysis-outputs", "data"),
        State("store-databases-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_annotated_tree(phylo_data, analysis_data, databases_data, run_config_data):
        from phylo.phylo_diversity import load_asv_mapping, tree_layout
        from phylo.phylo_metrics import support_tree_path
        from ui.plots import annotated_tree_figure

        from ui.feature_sets import phylo_dir, resolve_feature_set

        files = _phylo_files(phylo_data)
        output_dir = (run_config_data or {}).get("output_dir")
        # Colour nodes by support, so read the tree that actually carries it — from the
        # per-set phylo dir, since the on-disk support file lives under phylo/<set>/raxml/.
        feature_set = resolve_feature_set(None, _normalize_store(analysis_data).get("outputs") or {}, output_dir)
        tree_path = support_tree_path(files, phylo_dir(output_dir, feature_set) if output_dir else None)
        if tree_path is None:
            return html.Div("Run the phylo stage to see the tree.", className="empty-state-sub")

        layout = tree_layout(tree_path, max_tips=250)
        if layout is None:
            return html.Div("Tree too large to render in-app (>250 tips); download it from Results.", className="mini-note")

        mapping_csv = _mapping_csv(files, run_config_data)
        mapping = load_asv_mapping(mapping_csv) if mapping_csv else {}
        taxonomy = _tip_taxonomy(mapping, databases_data, run_config_data)
        abundance = _tip_abundance(mapping, analysis_data, run_config_data)

        figure = annotated_tree_figure(layout, taxonomy=taxonomy, abundance=abundance)
        if figure is None:
            return html.Div("Tree could not be rendered.", className="empty-state-sub")
        children = [dcc.Graph(figure=figure, config=_GRAPH_CONFIG)]
        if not mapping:
            children.append(
                html.Div(
                    "No ASV id mapping found, so tips could not be linked to the abundance table or taxonomy.",
                    className="mini-note",
                    style={"color": "var(--amber)"},
                )
            )
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "8px"})

    @app.callback(
        Output("phylo-support", "children"),
        Input("store-phylo-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_branch_support(phylo_data, run_config_data):
        import pandas as pd

        from ui.plots import branch_support_figure

        paths = _phylo_metric_paths(phylo_data, run_config_data)
        if "branch_support_csv" not in paths:
            return html.Div("This tree carries no branch support values.", className="empty-state-sub")
        try:
            support = pd.read_csv(paths["branch_support_csv"])
        except Exception as exc:
            return html.Div(f"Branch support unavailable: {exc}", className="empty-state-sub")
        if support.empty:
            return html.Div("This tree carries no branch support values.", className="empty-state-sub")

        supported = int(support["is_supported"].sum())
        children = [
            html.Div(
                [
                    html.Div(f"{supported} of {len(support)} internal nodes supported", className="status-pill"),
                    html.Div(f"median support {support['support'].median():.1f}%", className="status-pill"),
                ],
                className="chip-list",
            )
        ]
        figure = branch_support_figure(support)
        if figure is not None:
            children.append(dcc.Graph(figure=figure, config=_GRAPH_CONFIG))
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "8px"})

    @app.callback(
        Output("phylo-unifrac", "children"),
        Input("store-phylo-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_unifrac(phylo_data, run_config_data):
        import pandas as pd

        from ui.plots import unifrac_heatmap_figure

        paths = _phylo_metric_paths(phylo_data, run_config_data)
        children = []
        for key, weighted in (("weighted_unifrac_csv", True), ("unweighted_unifrac_csv", False)):
            if key not in paths:
                continue
            try:
                figure = unifrac_heatmap_figure(pd.read_csv(paths[key]), weighted=weighted)
            except Exception:
                figure = None
            if figure is not None:
                children.append(dcc.Graph(figure=figure, config=_GRAPH_CONFIG))
        if not children:
            return html.Div(
                "UniFrac needs the tree, the abundance table and the ASV id mapping. Run the analysis and phylo "
                "stages with at least two samples.",
                className="empty-state-sub",
            )
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "12px"})

    @app.callback(
        Output("phylo-visuals", "children"),
        Input("store-phylo-outputs", "data"),
        Input("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_phylo_visuals(phylo_data, analysis_data, run_config_data):
        from phylo.phylo_diversity import faith_pd
        from ui.callbacks_results import _find_abundance_matrix
        from ui.plots import pd_figure

        files = _phylo_files(phylo_data)
        tree_path = files.get("tree_nwk") or files.get("best_tree") or files.get("support_tree")
        if not tree_path or not Path(tree_path).exists():
            return html.Div("Run the phylo stage to see phylogenetic diversity.", className="empty-state-sub")

        analysis = _normalize_store(analysis_data)
        matrix_path = _find_abundance_matrix(analysis.get("outputs") or {}, (run_config_data or {}).get("output_dir"))
        if not matrix_path:
            return html.Div("No abundance matrix: run the analysis stage.", className="empty-state-sub")

        faith = faith_pd(tree_path, matrix_path, id_mapping_csv=_mapping_csv(files, run_config_data))
        figure = pd_figure(faith)
        if figure is None:
            return html.Div("Phylogenetic diversity unavailable.", className="empty-state-sub")
        note = html.Div(
            f"{faith.get('matched_features', 0)} of {faith.get('total_tips', 0)} tree tips matched the abundance "
            "table. Faith's PD and richness measure different things: a sample can be rich in closely related "
            "features and still have little phylogenetic diversity.",
            className="mini-note",
        )
        return html.Div([note, dcc.Graph(figure=figure, config=_GRAPH_CONFIG)])

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
        seen_paths: set[str] = set()
        for key, label in [("tree_nwk", "Tree (.nwk)"), ("best_tree", "Best Tree (.nwk)"), ("support_tree", "Support Tree"), ("bootstraps", "Bootstraps"), ("best_model", "Best Model")]:
            path_str = files.get(key)
            if not path_str or path_str in seen_paths:
                continue
            p = Path(path_str)
            if not p.exists():
                continue
            seen_paths.add(path_str)
            try:
                rel = p.relative_to(runs_root)
                href = f"/runs/{rel.as_posix()}"
            except ValueError:
                href = "#"
            items.append(html.Div(html.A(f"⬇ {label}", href=href, target="_blank", className="btn-parafish")))

        if not items:
            return html.Div("No phylogenetic outputs available yet.", className="empty-state-sub")
        return html.Div(items, style={"display": "flex", "flexDirection": "column", "gap": "8px"})

def register_phylo_kpis(app) -> None:
    """Layer-1 KPI cards for Phylo, matching the other pages' stat-card contract."""
    from dash import Input, Output, State, html

    from ui.components.cards import status_badge

    @app.callback(
        Output("phylo-stat-tips", "children"),
        Output("phylo-stat-alignment", "children"),
        Output("phylo-stat-support", "children"),
        Input("store-phylo-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _phylo_kpis(_phylo, run_config_data, project_data):
        from pathlib import Path

        import pandas as pd

        from core.run_registry import load_run_state, run_dir_for, stage_record

        def card(value, label, hint=None, color=None):
            body = [
                html.Div(str(value), className="stat-value", style={"color": color} if color else {}),
                html.Div(label, className="stat-label"),
            ]
            if hint:
                body.append(html.Div(hint, className="mini-note"))
            return body

        run_dir = run_dir_for(run_config_data, project_data)
        state = load_run_state(run_dir)
        record = stage_record(state, "phylo")
        status = str(record.get("status") or "pending")
        if status in ("pending", "skipped") or not record.get("outputs"):
            blank = card("—", "Phylo", hint="a etapa não correu")
            return blank, blank, blank

        from ui.feature_sets import phylo_tree_path

        tips = "—"
        tree = phylo_tree_path(run_dir, "otu")
        if tree:
            tips = Path(tree).read_text(errors="replace").count(",") + 1

        metrics = {k: v for k, v in (record.get("metrics_files") or {}).items() if isinstance(v, str) and Path(v).exists()}

        align_len, support_median = "—", "—"
        alignment = metrics.get("otu.alignment_csv") or metrics.get("alignment_csv")
        if alignment:
            try:
                align_len = f"{len(pd.read_csv(alignment)):,} col"
            except Exception:
                pass
        support = metrics.get("otu.support_csv") or metrics.get("support_csv")
        if support:
            try:
                frame = pd.read_csv(support)
                if "support" in frame.columns and not frame.empty:
                    support_median = f"{frame['support'].median():.0f}"
            except Exception:
                pass

        return (
            card(tips, "Sequências na árvore"),
            card(align_len, "Comprimento do alinhamento"),
            card(support_median, "Suporte bootstrap (mediana)"),
        )
