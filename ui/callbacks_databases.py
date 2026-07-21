"""Callbacks related to the databases stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, dash_table, dcc, html, no_update
from dash.exceptions import PreventUpdate

from core.config import get_db_paths
from databases.bold_db import BOLDConfig
from databases.custom_db import CustomDBConfig
from databases.ncbi_db import NCBIConfig
from databases.pr2_db import PR2Config
from databases.reference_filter import resolve_marker
from databases.silva_db import SILVAConfig
from databases.taxonomy_runner import run_taxonomy_pipeline
from preprocess.cutadapt_runner import resolve_primers
from ui.components.collapsible import details
from ui.components.table import data_table
from ui.pipeline_state import _normalize_store, get_stage_record, stage_ready, stage_status, utc_now_iso
from ui.runtime_state import decode_upload_content, latest_run_dir, resolve_project_data, run_config_from_store, serialize_step_output


# Plain-language headers: "Best_Identity" is a column name, "% identidade" is what the
# reader is actually looking for.
_CONSENSUS_HEADERS = {
    "ASV_ID": "ID",
    "feature_type": "Tipo",
    "Best_Identity": "% identidade",
    "Best_Coverage": "% cobertura",
    "Best_Source": "Fonte",
    "Reliability": "Fiabilidade",
    "Consensus_Kingdom": "Reino",
    "Consensus_Phylum": "Filo",
    "Consensus_Class": "Classe",
    "Consensus_Order": "Ordem",
    "Consensus_Family": "Família",
    "Consensus_Genus": "Género",
    "Consensus_Species": "Espécie",
}



_ID_COLUMNS = ("ASV_ID", "asv_id", "qseqid", "id", "seq_id", "query", "feature")


def _clean_feature_ids(frame):
    """Strip VSEARCH/Illumina cruft from whatever column holds the feature id.

    VSEARCH's ``--relabel_keep`` writes ">OTU1;size=132828 A01426:1215:HW7YLDRX5:..."
    so the centroid read stays traceable — useful, but it means the raw source tables
    show a query id like ``OTU1;size=132828 A01426:...``. "OTU1" is the citable
    identifier; the rest is provenance, and it belongs in otu_id_map.csv, not in a
    column the user is trying to read.
    """
    for column in _ID_COLUMNS:
        if column in frame.columns:
            frame = frame.copy()
            frame[column] = (
                frame[column].astype(str).str.split(";").str[0].str.split().str[0]
            )
            break
    return frame


def _load_hits_frame(db_outputs: dict, run_output_dir, feature_set: str, min_identity: float, min_coverage: float, top_n: int):
    """The ranked hits table, filtered. Returns ``(frame, error)``."""
    import pandas as pd

    from ui.feature_sets import taxonomy_file

    path = taxonomy_file(db_outputs, run_output_dir, feature_set, "taxonomy_hits_csv")
    if not path:
        return None, None
    try:
        frame = pd.read_csv(path)
    except Exception as exc:
        return None, str(exc)
    if frame.empty:
        return frame, None

    for column in ("identity", "coverage", "hit_rank"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

    if min_identity and "identity" in frame.columns:
        frame = frame[frame["identity"].fillna(-1) >= float(min_identity)]
    if min_coverage and "coverage" in frame.columns:
        # A hit with no coverage (SINTAX) is kept when no coverage floor is asked for,
        # and dropped when one is — it cannot demonstrate it meets the bar.
        frame = frame[frame["coverage"].fillna(-1) >= float(min_coverage)]
    if top_n and "hit_rank" in frame.columns:
        frame = frame[frame["hit_rank"] <= int(top_n)]

    sort_columns = [c for c in ("ASV_ID", "hit_rank") if c in frame.columns]
    if sort_columns:
        frame = frame.sort_values(by=sort_columns)
    return frame, None


def _is_supported_custom_db_filename(filename: str) -> bool:
    name = (filename or "").lower()
    return name.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz"))


def _find_consensus_csv(db_outputs: dict, run_output_dir, feature_set: str = "otu") -> str | None:
    """The consensus table for ONE feature set. Never a blend of both."""
    from ui.feature_sets import taxonomy_file

    resolved = taxonomy_file(db_outputs, run_output_dir, feature_set, "taxonomy_consensus_csv")
    if resolved:
        return resolved
    # Pre-feature-set runs kept a single flat table.
    files = ((db_outputs or {}).get("taxonomy") or {}).get("files") or {}
    value = files.get("taxonomy_consensus_csv")
    if isinstance(value, str) and value and Path(value).exists():
        return value
    return None


def _find_source_status(db_outputs: dict, run_output_dir) -> str | None:
    files = ((db_outputs or {}).get("taxonomy") or {}).get("files") or {}
    value = files.get("source_status_csv")
    if isinstance(value, str) and value and Path(value).exists():
        return value
    if run_output_dir:
        from databases.database_metrics import taxonomy_table_path

        candidate = taxonomy_table_path(Path(run_output_dir) / "databases", "source_status.csv")
        if candidate:
            return str(candidate)
    return None


def _source_status_panel(status_csv: str):
    """Colour-coded per-source status so a crashed DB is never a silent 'None'."""
    import csv

    rows = []
    try:
        with open(status_csv, "r", encoding="utf-8", errors="replace") as handle:
            rows = list(csv.DictReader(handle))
    except Exception:
        return None
    if not rows:
        return None
    palette = {
        "ok": ("#3fb950", "rgba(63,185,80,0.12)", "✓ classificou"),
        "no_hits": ("#d29922", "rgba(210,153,34,0.12)", "sem atribuições"),
        # Neutral, not amber: a source that does not index this marker has not
        # under-performed, and colouring it like "found nothing" invites the reader to
        # conclude the taxa are absent from a database that was never asked.
        "not_applicable": ("#8b98aa", "rgba(139,152,170,0.12)", "marcador incompatível"),
        "skipped": ("#8b98aa", "rgba(139,152,170,0.12)", "ignorada"),
        "failed": ("#f85149", "rgba(248,81,73,0.12)", "falhou"),
    }
    pills = []
    for row in rows:
        status = (row.get("status") or "").strip().lower()
        color, bg, label = palette.get(status, ("#8b98aa", "rgba(139,152,170,0.12)", status or "—"))
        pills.append(
            html.Div(
                [
                    html.Div(row.get("source", ""), style={"fontWeight": "600", "color": "#f4f7fb"}),
                    html.Div(label, style={"color": color, "fontSize": "12px", "fontWeight": "600"}),
                    html.Div(row.get("detail", ""), className="mini-note"),
                ],
                style={"padding": "8px 10px", "borderRadius": "8px", "background": bg, "border": f"1px solid {color}55", "minWidth": "160px", "flex": "1"},
            )
        )
    return html.Div(
        [
            html.Div("Reference databases — status of this run", className="mini-note"),
            html.Div(pills, style={"display": "flex", "flexWrap": "wrap", "gap": "8px"}),
        ]
    )


def _find_custom_besthit(db_outputs: dict, run_output_dir) -> str | None:
    files = ((db_outputs or {}).get("taxonomy") or {}).get("files") or {}
    value = files.get("custom_blast_besthit_csv")
    if isinstance(value, str) and value and Path(value).exists():
        return value
    if run_output_dir:
        candidate = Path(run_output_dir) / "databases" / "custom_blast" / "custom_blast_besthit.csv"
        if candidate.exists():
            return str(candidate)
    return None



_GRAPH_CONFIG = {"displaylogo": False, "responsive": True}


def run_db_classification(*args, **kwargs):
    """Build a compact database summary from shared store data."""
    project_data = resolve_project_data(kwargs.get("project_data") or {})
    samples_data = kwargs.get("samples_data") or {}
    step_outputs = kwargs.get("step_outputs") or {}
    stage_record = get_stage_record(step_outputs, "databases")
    sample_count = len(samples_data.get("samples", []))
    project_name = project_data.get("project_name") or "No project selected"
    status = stage_status(stage_record)
    return html.Div(
        [
            html.Div("Database routing", className="status-pill"),
            html.Div(f"Stage status: {status}", className="status-pill"),
            html.Div(f"Project: {project_name}", className="mini-note"),
            html.Div(f"{sample_count} sample group(s) awaiting taxonomy consensus", className="mini-note"),
        ],
        className="chip-list",
    )


def _source_tables(db_outputs: dict) -> list[tuple[str, str]]:
    """Return ``(label, path)`` for each per-source table the taxonomy stage wrote."""
    files = ((db_outputs or {}).get("taxonomy") or {}).get("files") or {}
    labels = {
        "silva_sintax_csv": "SILVA (SINTAX)",
        "pr2_sintax_csv": "PR2 (SINTAX)",
        "custom_sintax_csv": "Custom DB (SINTAX)",
        "custom_blast_besthit_csv": "Custom DB (BLAST best hit)",
        "blast_tsv": "NCBI (BLAST)",
        "bold_xlsx": "BOLD",
        "conflict_report_csv": "Conflicts between sources",
    }
    found = []
    for key, label in labels.items():
        value = files.get(key)
        if isinstance(value, str) and value and Path(value).exists() and Path(value).stat().st_size > 0:
            found.append((label, value))
    return found


def _read_preview(path: str, max_rows: int = 200):
    """Read a delimited or spreadsheet table into (columns, rows) for preview."""
    import pandas as pd

    suffix = Path(path).suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        frame = pd.read_excel(path, nrows=max_rows)
    else:
        sep = "\t" if suffix in {".tsv", ".txt"} else ","
        # VSEARCH's SINTAX --tabbedout is a headerless TSV written to a .csv name.
        if "sintax" in Path(path).name:
            frame = pd.read_csv(path, sep="\t", header=None, nrows=max_rows)
            frame.columns = ["Query", "Taxonomy (bootstrap)", "Strand", "Taxonomy (above cutoff)"][: len(frame.columns)]
        else:
            frame = pd.read_csv(path, sep=sep, nrows=max_rows)
    return frame


def _resolve_output_dir(run_config_data, project_data=None) -> str | None:
    """The run's output directory, from the session store or reconstructed from disk.

    ``run_config_data`` is per-session: a fresh tab or a restarted app has it empty, so
    ``output_dir`` is gone even though the run is still on disk. Fall back to the project's
    newest run (``latest_run_dir``, which now orders by run_id, not mtime) instead of
    blanking the page — the disk path the Results page already takes. This is Eixo B only:
    it resolves WHICH run, never which feature set.
    """
    output_dir = (run_config_data or {}).get("output_dir")
    if output_dir:
        return str(output_dir)
    run_dir = latest_run_dir(project_data) if project_data else None
    return str(run_dir) if run_dir else None


def _database_metric_paths(db_data, run_config_data, project_data=None) -> dict[str, str]:
    """Locate the taxonomy tables the runner persisted."""
    record = _normalize_store(db_data)
    recorded = record.get("metrics_files") or {}
    paths = {key: value for key, value in recorded.items() if isinstance(value, str) and Path(value).exists()}
    if paths:
        return paths
    output_dir = _resolve_output_dir(run_config_data, project_data)
    if not output_dir:
        return {}
    from databases.database_metrics import database_metrics_paths

    return {key: str(path) for key, path in database_metrics_paths(Path(output_dir) / "databases").items()}


def _taxonomy_file(db_outputs: dict, run_output_dir, key: str, fallback: str) -> str | None:
    """Path to one taxonomy artifact, from the record or its conventional location."""
    files = ((db_outputs or {}).get("taxonomy") or {}).get("files") or {}
    value = files.get(key)
    if isinstance(value, str) and value and Path(value).exists():
        return value
    if run_output_dir:
        from databases.database_metrics import taxonomy_table_path

        candidate = taxonomy_table_path(Path(run_output_dir) / "databases", fallback)
        if candidate:
            return str(candidate)
    return None


def register_callbacks(app: Any) -> None:
    """Register database callbacks in the Dash application."""

    @app.callback(
        Output("databases-feature-set", "options"),
        Output("databases-feature-set", "value"),
        Output("databases-feature-set-note", "children"),
        Input("store-analysis-outputs", "data"),
        State("databases-feature-set", "value"),
        State("store-run-config", "data"),
    )
    def _sync_feature_sets(analysis_out, current, run_config_data):
        """Offer only the feature sets this run actually produced."""
        from ui.feature_sets import available_feature_sets

        analysis = _normalize_store(analysis_out)
        outputs = analysis.get("outputs") or {}
        run_output_dir = (run_config_data or {}).get("output_dir")
        sets = available_feature_sets(outputs, run_output_dir)
        options = [{"label": f" {entry['label']}", "value": entry["key"]} for entry in sets]
        keys = [entry["key"] for entry in sets]
        value = current if current in keys else (keys[0] if keys else "otu")
        if len(keys) > 1:
            note = "Este run gerou os dois conjuntos. Alterna aqui; para os comparar, vê o painel «OTU vs ASV»."
        elif keys:
            note = f"Este run gerou apenas {sets[0]['label']}. Escolhe o método «both» no Setup para gerar ambos."
        else:
            note = ""
        return options, value, note

    @app.callback(
        Output("databases-hits", "children"),
        Input("store-databases-outputs", "data"),
        Input("databases-feature-set", "value"),
        Input("hits-min-identity", "value"),
        Input("hits-min-coverage", "value"),
        Input("hits-top-n", "value"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_hits(db_out, feature_set, min_identity, min_coverage, top_n, run_config_data, project_data):
        db = _normalize_store(db_out)
        outputs = db.get("outputs") or {}
        run_output_dir = _resolve_output_dir(run_config_data, project_data)
        feature_set = feature_set or "otu"

        frame, error = _load_hits_frame(outputs, run_output_dir, feature_set, min_identity or 0, min_coverage or 0, top_n or 0)
        if error:
            return html.Div(f"Tabela de hits indisponível: {error}", className="empty-state-sub")
        if frame is None:
            return html.Div(
                "Ainda não há hits. Corre a etapa Databases com pelo menos uma fonte que alinhe "
                "(NCBI, BOLD ou base customizada) — o SINTAX não produz % identidade.",
                className="empty-state-sub",
            )
        if frame.empty:
            return html.Div("Nenhum hit passa os filtros escolhidos.", className="empty-state-sub")

        n_features = frame["ASV_ID"].nunique() if "ASV_ID" in frame.columns else 0
        headers = {**_CONSENSUS_HEADERS, "hit_rank": "Hit #", "identity": "% identidade", "coverage": "% cobertura", "reliability": "Fiabilidade", "source": "Fonte", "taxon": "Taxon", "subject": "Referência"}
        return html.Div(
            [
                html.Div(
                    [
                        html.Div(f"{len(frame)} hits", className="status-pill"),
                        html.Div(f"{n_features} features", className="status-pill"),
                    ],
                    className="chip-list",
                ),
                data_table(
                    data=frame.to_dict("records"),
                    columns=[{"name": headers.get(str(c), str(c)), "id": str(c)} for c in frame.columns],
                    sort_action="native",
                    filter_action="native",
                    row_selectable="multi",
                    id="hits-table",
                    page_size=15,
                ),
            ],
            style={"display": "flex", "flexDirection": "column", "gap": "10px"},
        )

    @app.callback(
        Output("download-hits-csv", "data"),
        Input("btn-download-hits", "n_clicks"),
        State("store-databases-outputs", "data"),
        State("databases-feature-set", "value"),
        State("hits-min-identity", "value"),
        State("hits-min-coverage", "value"),
        State("hits-top-n", "value"),
        State("store-run-config", "data"),
        State("store-project", "data"),
        prevent_initial_call=True,
    )
    def _download_hits(n_clicks, db_out, feature_set, min_identity, min_coverage, top_n, run_config_data, project_data):
        """Export exactly what the user is looking at, as one flat table."""
        if not n_clicks:
            raise PreventUpdate
        db = _normalize_store(db_out)
        outputs = db.get("outputs") or {}
        run_output_dir = _resolve_output_dir(run_config_data, project_data)
        feature_set = feature_set or "otu"
        frame, error = _load_hits_frame(outputs, run_output_dir, feature_set, min_identity or 0, min_coverage or 0, top_n or 0)
        if frame is None or error or frame.empty:
            raise PreventUpdate
        return dcc.send_data_frame(frame.to_csv, f"{feature_set}_best_hits.csv", index=False)

    @app.callback(
        Output("databases-feature-comparison", "children"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_feature_comparison(db_out, run_config_data, project_data):
        import pandas as pd

        db = _normalize_store(db_out)
        outputs = db.get("outputs") or {}
        comparison = outputs.get("feature_set_comparison") or {}
        run_output_dir = _resolve_output_dir(run_config_data, project_data)

        summary_path = comparison.get("feature_set_summary_csv")
        taxa_path = comparison.get("feature_set_taxa_csv")
        if not summary_path and run_output_dir:
            candidate = Path(run_output_dir) / "analysis" / "comparison" / "feature_set_summary.csv"
            summary_path = str(candidate) if candidate.exists() else None
            candidate = Path(run_output_dir) / "analysis" / "comparison" / "feature_set_taxa.csv"
            taxa_path = str(candidate) if candidate.exists() else None

        if not summary_path or not Path(summary_path).exists():
            return html.Div(
                "Disponível quando o run gera os dois conjuntos. Escolhe o método de análise "
                "«both» (VSEARCH + DADA2) no Setup.",
                className="empty-state-sub",
            )

        children: list = []
        try:
            summary = pd.read_csv(summary_path)
        except Exception as exc:
            return html.Div(f"Comparação indisponível: {exc}", className="empty-state-sub")

        headers = {
            "feature_set": "Conjunto",
            "label": "Método",
            "n_features": "Nº features",
            "n_samples": "Nº amostras",
            "total_reads": "Reads",
            "features_in_multiple_samples": "Em várias amostras",
            "features_in_one_sample": "Só numa amostra",
            "mean_reads_per_feature": "Reads/feature",
        }
        children.append(html.Div("Resumo — uma linha por conjunto, nunca uma linha com os dois", className="mini-note"))
        children.append(
            data_table(
                data=summary.to_dict("records"),
                columns=[{"name": headers.get(str(c), str(c)), "id": str(c)} for c in summary.columns],
            )
        )

        if taxa_path and Path(taxa_path).exists():
            try:
                taxa = pd.read_csv(taxa_path)
            except Exception:
                taxa = None
            if taxa is not None and not taxa.empty:
                both = int((taxa["agreement"] == "both").sum()) if "agreement" in taxa.columns else 0
                otu_only = int((taxa["agreement"] == "otu_only").sum()) if "agreement" in taxa.columns else 0
                asv_only = int((taxa["agreement"] == "asv_only").sum()) if "agreement" in taxa.columns else 0
                taxa_headers = {
                    "taxon": "Taxon",
                    "rank": "Nível",
                    "n_otus": "Nº OTUs",
                    "n_asvs": "Nº ASVs",
                    "otu_reads": "Reads (OTU)",
                    "asv_reads": "Reads (ASV)",
                    "agreement": "Concordância",
                }
                children.append(
                    html.Div(
                        [
                            html.Div(f"{both} taxa detetados pelos dois", className="status-pill"),
                            html.Div(f"{otu_only} só por OTUs", className="status-pill"),
                            html.Div(f"{asv_only} só por ASVs", className="status-pill"),
                        ],
                        className="chip-list",
                        style={"marginTop": "12px"},
                    )
                )
                children.append(
                    html.Div(
                        "Ponte por taxon: cada lado mantém a sua própria contagem de features e de reads, "
                        "em colunas separadas. Nunca são somadas.",
                        className="mini-note",
                    )
                )
                children.append(
                    data_table(
                        data=taxa.to_dict("records"),
                        columns=[
                            {"name": taxa_headers.get(str(c), str(c)), "id": str(c)}
                            for c in taxa.columns
                            if c not in ("in_otu", "in_asv")
                        ],
                        sort_action="native",
                        filter_action="native",
                        page_size=12,
                        conditional=[
                            {"if": {"filter_query": '{agreement} = "both"', "column_id": "agreement"}, "color": "#3fb950"},
                            {"if": {"filter_query": '{agreement} = "otu_only"', "column_id": "agreement"}, "color": "#d29922"},
                            {"if": {"filter_query": '{agreement} = "asv_only"', "column_id": "agreement"}, "color": "#58a6ff"},
                        ],
                    )
                )

        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "10px"})

    @app.callback(
        Output("databases-resolution", "children"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_resolution(db_out, run_config_data, project_data):
        import pandas as pd

        from ui.plots import taxonomic_resolution_figure, taxonomy_sunburst_figure

        paths = _database_metric_paths(db_out, run_config_data, project_data)
        db = _normalize_store(db_out)
        run_output_dir = _resolve_output_dir(run_config_data, project_data)
        children = []

        if "resolution_csv" in paths:
            try:
                resolution = pd.read_csv(paths["resolution_csv"])
            except Exception as exc:
                return html.Div(f"Resolution unavailable: {exc}", className="empty-state-sub")
            figure = taxonomic_resolution_figure(resolution)
            if figure is not None:
                children.append(dcc.Graph(figure=figure, config=_GRAPH_CONFIG))
            children.append(
                data_table(
                    data=resolution.to_dict("records"),
                    columns=[
                        {"name": "Rank", "id": "rank"},
                        {"name": "Features", "id": "n_features", "type": "numeric"},
                        {"name": "% features", "id": "pct_features"},
                        {"name": "Reads", "id": "n_reads", "type": "numeric", "format": {"specifier": ","}},
                        {"name": "% reads", "id": "pct_reads"},
                    ],
                    page_size=9,
                )
            )

        consensus_path = _find_consensus_csv(db.get("outputs") or {}, run_output_dir)
        if consensus_path:
            try:
                sunburst = taxonomy_sunburst_figure(pd.read_csv(consensus_path))
            except Exception:
                sunburst = None
            if sunburst is not None:
                children.append(dcc.Graph(figure=sunburst, config=_GRAPH_CONFIG))

        if not children:
            return html.Div("Run the Databases stage to see taxonomic resolution.", className="empty-state-sub")
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "12px"})

    @app.callback(
        Output("databases-confidence", "children"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_confidence(db_out, run_config_data, project_data):
        import pandas as pd

        from ui.plots import confidence_histogram_figure, reliability_figure

        db = _normalize_store(db_out)
        run_output_dir = _resolve_output_dir(run_config_data, project_data)
        assignments_path = _taxonomy_file(db.get("outputs") or {}, run_output_dir, "source_assignments_csv", "source_assignments.csv")
        if not assignments_path:
            return html.Div("Run the Databases stage to see per-assignment confidence.", className="empty-state-sub")
        try:
            assignments = pd.read_csv(assignments_path)
        except Exception as exc:
            return html.Div(f"Assignments unavailable: {exc}", className="empty-state-sub")
        if assignments.empty:
            return html.Div("No taxonomic assignments were made.", className="empty-state-sub")

        children = []
        for figure in (confidence_histogram_figure(assignments), reliability_figure(assignments)):
            if figure is not None:
                children.append(dcc.Graph(figure=figure, config=_GRAPH_CONFIG))

        red, amber, green = "#f85149", "#d29922", "#3fb950"
        children.append(
            data_table(
                data=assignments.to_dict("records"),
                columns=[
                    {"name": "OTU / ASV", "id": "ASV_ID"},
                    {"name": "Source", "id": "source"},
                    {"name": "Deepest rank", "id": "rank"},
                    {"name": "Taxon", "id": "taxon"},
                    {"name": "% identity", "id": "identity"},
                    {"name": "% coverage", "id": "coverage"},
                    {"name": "Confidence", "id": "confidence"},
                    {"name": "Reliability", "id": "reliability"},
                ],
                sort_action="native",
                filter_action="native",
                page_size=12,
                conditional=[
                    {"if": {"filter_query": '{reliability} = "low"', "column_id": "reliability"}, "color": red, "fontWeight": "600"},
                    {"if": {"filter_query": '{reliability} = "moderate"', "column_id": "reliability"}, "color": amber, "fontWeight": "600"},
                    {"if": {"filter_query": '{reliability} = "high"', "column_id": "reliability"}, "color": green, "fontWeight": "600"},
                ],
            )
        )
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "12px"})

    @app.callback(
        Output("databases-conflicts", "children"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_conflicts(db_out, run_config_data, project_data):
        import pandas as pd

        db = _normalize_store(db_out)
        run_output_dir = _resolve_output_dir(run_config_data, project_data)
        path = _taxonomy_file(db.get("outputs") or {}, run_output_dir, "conflict_report_csv", "conflict_report.csv")
        if not path:
            return html.Div("Run the Databases stage to compare sources.", className="empty-state-sub")
        try:
            conflicts = pd.read_csv(path)
        except Exception as exc:
            return html.Div(f"Conflict report unavailable: {exc}", className="empty-state-sub")
        if conflicts.empty:
            return html.Div(
                "No conflicts: every OTU/ASV classified by more than one source got the same answer at every rank.",
                className="mini-note",
                style={"color": "var(--green)"},
            )
        return html.Div(
            [
                html.Div(
                    f"{len(conflicts)} OTU/ASV(s) where the sources disagree. The consensus falls back to the "
                    "deepest rank they do agree on.",
                    className="mini-note",
                ),
                data_table(
                    data=conflicts.to_dict("records"),
                    columns=[{"name": str(column), "id": str(column)} for column in conflicts.columns],
                    sort_action="native",
                    filter_action="native",
                    page_size=12,
                ),
            ]
        )

    @app.callback(
        Output("databases-contaminants", "children"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_contaminants(db_out, run_config_data, project_data):
        import pandas as pd

        from ui.plots import contaminant_figure

        paths = _database_metric_paths(db_out, run_config_data, project_data)
        if "contaminants_csv" not in paths:
            return html.Div(
                "No negative controls were marked. Tick 'Control' on the blank samples in the Samples table to "
                "screen for reagent and cross contamination.",
                className="empty-state-sub",
            )
        try:
            contaminants = pd.read_csv(paths["contaminants_csv"])
        except Exception as exc:
            return html.Div(f"Contamination screen unavailable: {exc}", className="empty-state-sub")

        n_flagged = int(contaminants["is_contaminant"].sum())
        children = [
            html.Div(
                [
                    html.Div(f"{n_flagged} feature(s) flagged as contaminants", className="status-pill"),
                    html.Div(f"{len(contaminants)} feature(s) screened", className="status-pill"),
                ],
                className="chip-list",
            ),
            html.Div(
                "Flagged by prevalence, not abundance: a feature that turns up in the blanks and few real samples "
                "is contamination however abundant it is. With only a couple of controls the test has little power, "
                "so an unflagged feature is not evidence of cleanliness.",
                className="mini-note",
            ),
        ]
        figure = contaminant_figure(contaminants)
        if figure is not None:
            children.append(dcc.Graph(figure=figure, config=_GRAPH_CONFIG))
        children.append(
            data_table(
                data=contaminants.head(200).to_dict("records"),
                columns=[
                    {"name": "Feature", "id": "feature"},
                    {"name": "Prev. controls", "id": "prevalence_controls"},
                    {"name": "Prev. samples", "id": "prevalence_samples"},
                    {"name": "Reads (controls)", "id": "reads_controls", "type": "numeric"},
                    {"name": "Reads (samples)", "id": "reads_samples", "type": "numeric"},
                    {"name": "Score", "id": "score"},
                    {"name": "Contaminant", "id": "is_contaminant"},
                ],
                sort_action="native",
                page_size=12,
                conditional=[
                    {"if": {"filter_query": "{is_contaminant} = 1"}, "backgroundColor": "rgba(248,81,73,0.10)"},
                ],
            )
        )
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "10px"})

    @app.callback(
        Output("databases-prevalence", "children"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_prevalence(db_out, run_config_data, project_data):
        import pandas as pd

        from ui.plots import prevalence_figure

        paths = _database_metric_paths(db_out, run_config_data, project_data)
        if "taxon_prevalence_csv" not in paths:
            return html.Div("Run the Databases stage to see per-taxon prevalence.", className="empty-state-sub")
        try:
            prevalence = pd.read_csv(paths["taxon_prevalence_csv"])
        except Exception as exc:
            return html.Div(f"Prevalence unavailable: {exc}", className="empty-state-sub")

        children = []
        if "permanova_csv" in paths:
            try:
                result = pd.read_csv(paths["permanova_csv"]).iloc[0]
                children.append(
                    html.Div(
                        f"{result['test']}: pseudo-F = {result['pseudo_f']}, p = {result['p_value']} "
                        f"({int(result['n_groups'])} groups, {int(result['n_samples'])} samples, "
                        f"{int(result['permutations'])} permutations).",
                        className="status-pill",
                    )
                )
            except Exception:
                pass
        figure = prevalence_figure(prevalence)
        if figure is not None:
            children.append(dcc.Graph(figure=figure, config=_GRAPH_CONFIG))
        children.append(
            data_table(
                data=prevalence.to_dict("records"),
                columns=[
                    {"name": "Taxon", "id": "taxon"},
                    {"name": "Detected in", "id": "n_detected", "type": "numeric"},
                    {"name": "Samples", "id": "n_samples", "type": "numeric"},
                    {"name": "Prevalence", "id": "prevalence"},
                    {"name": "95% CI low", "id": "ci_low"},
                    {"name": "95% CI high", "id": "ci_high"},
                    {"name": "Total reads", "id": "total_reads", "type": "numeric", "format": {"specifier": ","}},
                ],
                sort_action="native",
                filter_action="native",
                page_size=15,
            )
        )
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "10px"})

    @app.callback(
        Output("databases-summary", "children", allow_duplicate=True),
        Input("store-databases-outputs", "data"),
        State("store-project", "data"),
        State("store-samples", "data"),
        State("store-analysis-outputs", "data"),
        prevent_initial_call="initial_duplicate",
    )
    def _sync_databases_summary(db_out, project_data, samples_data, analysis_out):
        return run_db_classification(
            project_data=project_data,
            samples_data=samples_data,
            step_outputs={"analysis": _normalize_store(analysis_out), "databases": _normalize_store(db_out)},
        )

    @app.callback(
        Output("databases-consensus", "children"),
        Input("store-databases-outputs", "data"),
        Input("databases-feature-set", "value"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_consensus(db_out, feature_set, run_config_data, project_data):
        from ui.feature_sets import FEATURE_SET_NOUN

        db = _normalize_store(db_out)
        outputs = db.get("outputs") or {}
        run_output_dir = _resolve_output_dir(run_config_data, project_data)
        feature_set = feature_set or "otu"
        noun = FEATURE_SET_NOUN.get(feature_set, "feature")

        children: list = []
        consensus_path = _find_consensus_csv(outputs, run_output_dir, feature_set)
        if not consensus_path:
            children.append(
                html.Div(
                    f"Ainda não há consenso de taxonomia para {noun}s. Ativa pelo menos uma base de "
                    f"referência e corre a etapa Databases (e o método de análise que gera {noun}s).",
                    className="empty-state-sub",
                )
            )
        else:
            try:
                frame = _read_preview(consensus_path, max_rows=500)
            except Exception as exc:
                children.append(html.Div(f"Consensus table unavailable: {exc}", className="empty-state-sub"))
            else:
                classified = 0
                if "Consensus_Genus" in frame.columns:
                    classified = int(frame["Consensus_Genus"].astype(str).str.strip().replace("nan", "").astype(bool).sum())
                pills = [
                    html.Div(f"{len(frame)} {noun}s", className="status-pill"),
                    html.Div(f"{classified} classificados até género", className="status-pill"),
                ]
                if "Reliability" in frame.columns:
                    high = int((frame["Reliability"].astype(str) == "high").sum())
                    pills.append(html.Div(f"{high} de fiabilidade alta", className="status-pill"))
                children.append(html.Div(pills, className="chip-list"))

                # Sort by the evidence, best first, so the strongest calls are the ones
                # you see without scrolling.
                if "Best_Identity" in frame.columns:
                    frame = frame.sort_values(
                        by=[c for c in ("Best_Identity", "Best_Coverage") if c in frame.columns],
                        ascending=False,
                        na_position="last",
                    )
                children.append(
                    data_table(
                        data=frame.to_dict("records"),
                        columns=[{"name": _CONSENSUS_HEADERS.get(str(c), str(c)), "id": str(c)} for c in frame.columns],
                        sort_action="native",
                        filter_action="native",
                        page_size=12,
                    )
                )

        sources = _source_tables(outputs)
        if sources:
            # These are the RAW per-source dumps, straight from BLAST/SINTAX/BOLD. They
            # still carry unprocessed Illumina headers in the query column
            # ("OTU1;size=132828 A01426:1215:HW7YLDRX5:..."), which is an internal
            # artefact and not something to show a user by default. Tabs are right here
            # (the sources are mutually exclusive views of the same features), but the
            # whole block belongs behind a disclosure until the ids are clean.
            tabs = []
            for label, path in sources:
                try:
                    frame = _read_preview(path)
                except Exception as exc:
                    body = html.Div(f"Could not read {Path(path).name}: {exc}", className="empty-state-sub")
                else:
                    frame = _clean_feature_ids(frame)
                    body = data_table(
                        data=frame.to_dict("records"),
                        columns=[{"name": str(c), "id": str(c)} for c in frame.columns],
                        sort_action="native",
                        page_size=10,
                    )
                tabs.append(dcc.Tab(label=label, children=html.Div(body, style={"paddingTop": "10px"})))
            children.append(
                details(
                    "Resultados por fonte (dados em bruto)",
                    dcc.Tabs(tabs),
                    subtitle="saída directa de cada base — para auditoria",
                )
            )

        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "10px"})

    @app.callback(
        Output("databases-visuals", "children"),
        Input("store-databases-outputs", "data"),
        Input("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_db_visuals(db_out, analysis_out, run_config_data, project_data):
        from analysis.community_analysis import bray_curtis_matrix, pcoa, sample_summary, taxa_composition
        from ui.callbacks_results import _find_abundance_matrix
        from ui.plots import beta_heatmap_figure, composition_figure, ordination_figure

        db = _normalize_store(db_out)
        analysis = _normalize_store(analysis_out)
        run_output_dir = _resolve_output_dir(run_config_data, project_data)

        graph_config = {"displaylogo": False, "responsive": True}
        children: list = []

        # Per-source status first, so a crashed/skipped reference DB is explicit.
        status_csv = _find_source_status(db.get("outputs") or {}, run_output_dir)
        if status_csv:
            panel = _source_status_panel(status_csv)
            if panel is not None:
                children.append(panel)

        matrix_path = _find_abundance_matrix(analysis.get("outputs") or {}, run_output_dir)
        if not matrix_path:
            if children:
                children.append(html.Div("Run the analysis and databases stages to see taxonomic composition and sample comparison.", className="empty-state-sub"))
                return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "14px"})
            return html.Div("Run the analysis and databases stages to see taxonomic composition and sample comparison.", className="empty-state-sub")
        consensus_path = _find_consensus_csv(db.get("outputs") or {}, run_output_dir)

        comp_fig = composition_figure(taxa_composition(matrix_path, consensus_path, rank="genus", top_n=12))
        if comp_fig is not None:
            children.append(dcc.Graph(figure=comp_fig, config=graph_config))

        bc = bray_curtis_matrix(matrix_path)
        ord_fig = ordination_figure(pcoa(bc))
        if ord_fig is not None:
            children.append(dcc.Graph(figure=ord_fig, config=graph_config))
        beta_fig = beta_heatmap_figure(bc)
        if beta_fig is not None:
            children.append(dcc.Graph(figure=beta_fig, config=graph_config))

        summary = sample_summary(matrix_path, consensus_path, rank="genus")
        if summary:
            summary_sorted = sorted(summary, key=lambda s: s["reads"], reverse=True)
            children.append(
                html.Div(
                    [
                        html.Div("Per-sample summary — sortable to find the best/worst samples", className="mini-note"),
                        data_table(
                            data=summary_sorted,
                            columns=[{"name": c, "id": c} for c in ("sample", "reads", "richness", "shannon", "pct_classified")],
                            sort_action="native",
                            page_size=15,
                        ),
                    ]
                )
            )

        if not children:
            return html.Div("Not enough data for taxonomic comparison yet.", className="empty-state-sub")
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "14px"})

    @app.callback(
        Output("custom-comparison", "children"),
        Input("store-databases-outputs", "data"),
        Input("custom-cmp-minid", "value"),
        Input("custom-cmp-species", "value"),
        State("store-analysis-outputs", "data"),
        State("store-run-config", "data"),
    )
    def _render_custom_comparison(db_out, min_id, species_filter, analysis_out, run_config_data):
        from analysis.custom_db_compare import build_custom_comparison
        from ui.callbacks_results import _find_abundance_matrix
        from ui.plots import detection_heatmap_figure, match_rate_figure

        db = _normalize_store(db_out)
        analysis = _normalize_store(analysis_out)
        run_output_dir = (run_config_data or {}).get("output_dir")
        besthit = _find_custom_besthit(db.get("outputs") or {}, run_output_dir)
        if not besthit:
            return html.Div(
                "No custom-database comparison yet. Upload a reference FASTA above, then run the Databases stage.",
                className="empty-state-sub",
            )
        matrix_path = _find_abundance_matrix(analysis.get("outputs") or {}, run_output_dir)
        if not matrix_path:
            return html.Div("Run the analysis stage first to build the abundance matrix.", className="empty-state-sub")
        consensus_path = _find_consensus_csv(db.get("outputs") or {}, run_output_dir)

        comp = build_custom_comparison(
            besthit,
            matrix_path,
            min_identity=float(min_id or 0),
            species_filter=(species_filter or None),
            consensus_csv=consensus_path,
        )
        if comp.get("error"):
            return html.Div(comp["error"], className="empty-state-sub")

        summary = comp.get("summary") or {}
        graph_config = {"displaylogo": False, "responsive": True}
        children: list = []

        chips = html.Div(
            [
                html.Div(f"{summary.get('n_otus', 0)} OTUs/ASVs", className="status-pill"),
                html.Div(f"{summary.get('n_matched', 0)} matched", className="status-pill"),
                html.Div(f"{summary.get('n_unmatched', 0)} None (no match)", className="status-pill"),
                html.Div(f"{summary.get('pct_matched', 0)}% matched", className="status-pill"),
                html.Div(f"{summary.get('n_taxa', 0)} taxa · {summary.get('n_samples', 0)} samples", className="status-pill"),
            ],
            className="chip-list",
        )
        children.append(chips)

        donut = match_rate_figure(summary)
        if donut is not None:
            children.append(dcc.Graph(figure=donut, config=graph_config))
        heat = detection_heatmap_figure(comp.get("detection") or {})
        if heat is not None:
            children.append(dcc.Graph(figure=heat, config=graph_config))

        per_otu = comp.get("per_otu") or []
        max_rows = 1000
        truncated = len(per_otu) > max_rows
        table_rows = per_otu[:max_rows]
        columns = [
            {"name": "OTU / ASV", "id": "OTU_ID"},
            {"name": "Custom DB match", "id": "species"},
            {"name": "Public DBs (consensus)", "id": "public_species"},
            {"name": "% identity", "id": "pident"},
            {"name": "Total reads", "id": "total_reads"},
            {"name": "# samples", "id": "n_samples_detected"},
            {"name": "Status", "id": "status"},
        ]
        red = "#f85149"
        style_conditional = [
            {"if": {"filter_query": '{status} = "None"'}, "backgroundColor": "rgba(248,81,73,0.10)"},
            {"if": {"filter_query": '{species} = "None"', "column_id": "species"}, "color": red, "fontWeight": "600"},
            {"if": {"filter_query": '{public_species} = "None"', "column_id": "public_species"}, "color": red, "fontWeight": "600"},
            {"if": {"filter_query": '{status} = "None"', "column_id": "status"}, "color": red, "fontWeight": "600"},
            {"if": {"filter_query": '{status} = "Match"', "column_id": "status"}, "color": "#3fb950", "fontWeight": "600"},
        ]
        children.append(
            html.Div(
                [
                    html.Div(
                        "Per-OTU comparison — Custom DB match vs public-database consensus. 'None' (red) = no match. "
                        + (f"Showing top {max_rows} of {len(per_otu)} OTUs by reads." if truncated else "Sortable and filterable."),
                        className="mini-note",
                    ),
                    data_table(
                        data=table_rows,
                        columns=columns,
                        sort_action="native",
                        filter_action="native",
                        page_size=15,
                        conditional=style_conditional,
                    ),
                ]
            )
        )
        return html.Div(children, style={"display": "flex", "flexDirection": "column", "gap": "14px"})

    @app.callback(
        Output("custom-db-upload-status", "children"),
        Output("store-run-config", "data", allow_duplicate=True),
        Input("upload-custom-db", "contents"),
        State("upload-custom-db", "filename"),
        State("store-project", "data"),
        State("store-run-config", "data"),
        prevent_initial_call=True,
    )
    def _handle_custom_db_upload(contents, filename, project_data, run_config_data):
        if not contents or not filename:
            raise PreventUpdate
        if not _is_supported_custom_db_filename(filename):
            return html.Div(f"Unsupported format: {filename}. Use .fasta, .fa or .fasta.gz", style={"color": "red"}), no_update
        project = resolve_project_data(project_data)
        project_name = project.get("project_name") or "project"
        db_dir = Path("projects") / project_name / "custom_dbs"
        db_dir.mkdir(parents=True, exist_ok=True)
        dest = db_dir / Path(filename).name
        dest.write_bytes(decode_upload_content(contents))
        config = CustomDBConfig(dest)
        valid, error = config.validate()
        if not valid:
            dest.unlink(missing_ok=True)
            return html.Div(f"Ficheiro inválido: {error}", style={"color": "red"}), no_update
        raw = dict(run_config_data or {})
        params = raw.setdefault("vsearch_params", {})
        params["custom_db_path"] = str(dest)
        params["custom_db_kind"] = config.kind
        if config.kind == "blast":
            mode_note = "comparação por BLAST (blastn) contra os OTUs/ASVs"
        else:
            mode_note = "classificação SINTAX (cabeçalhos com ;tax=) integrada no consenso"
        return (
            html.Div(
                [
                    html.Div(f"✓ Carregada: {filename}", style={"color": "#3fb950", "fontWeight": "600"}),
                    html.Div(f"Modo detetado: {mode_note}.", className="mini-note"),
                    html.Div("Corre a etapa Databases para comparar as amostras contra esta base.", className="mini-note"),
                ]
            ),
            raw,
        )

    @app.callback(
        Output("databases-summary", "children"),
        Output("store-databases-outputs", "data"),
        Input("btn-run-databases", "n_clicks"),
        State("store-project", "data"),
        State("store-run-config", "data"),
        State("store-samples", "data"),
        State("store-analysis-outputs", "data"),
        State("store-databases-outputs", "data"),
        prevent_initial_call=True,
    )
    def _run_databases(n_clicks, project_data, run_config_data, samples_data, analysis_outputs, databases_outputs):
        if not n_clicks:
            raise PreventUpdate
        run_config = run_config_from_store(project_data, samples_data, run_config_data)
        analysis_record = _normalize_store(analysis_outputs)
        existing = _normalize_store(databases_outputs)
        signature = "|".join([run_config.run_id, run_config.branch.value, analysis_record.get("signature", "")])

        if existing.get("signature") == signature and stage_status(existing) in {"success", "failure", "blocked"}:
            return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": existing}), no_update

        if not stage_ready(analysis_record):
            record = {
                "signature": signature,
                "status": "blocked" if stage_status(analysis_record) == "failure" else "pending",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR" if stage_status(analysis_record) == "failure" else "INF", "msg": "Taxonomy aguardando análise concluída"}],
                "blocked_by": "analysis" if stage_status(analysis_record) == "failure" else None,
            }
            return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": record}), record

        combined_fasta_value = (analysis_record.get("outputs") or {}).get("combined_fasta")
        if not combined_fasta_value:
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": "Análise não produziu FASTA combinado"}],
                "error_message": "Análise não produziu FASTA combinado",
            }
            return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": record}), record

        fasta_path = Path(combined_fasta_value)
        if not fasta_path.exists():
            record = {
                "signature": signature,
                "status": "failure",
                "updated_at": utc_now_iso(),
                "outputs": {},
                "log_lines": [{"ts": utc_now_iso(), "level": "ERR", "msg": f"FASTA combinado nao encontrado: {fasta_path}"}],
                "error_message": f"FASTA combinado nao encontrado: {fasta_path}",
            }
            return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": record}), record
        db_paths = get_db_paths()
        custom_db_raw = (run_config.vsearch_params or {}).get("custom_db_path")
        custom_db_path = Path(custom_db_raw) if custom_db_raw and Path(custom_db_raw).exists() else None
        custom_db_kind = (run_config.vsearch_params or {}).get("custom_db_kind")
        from ui.callbacks_results import _find_abundance_matrix

        abundance_matrix_path = _find_abundance_matrix(analysis_record.get("outputs") or {}, run_config.output_dir)
        sintax_cutoff = float(run_config.vsearch_params.get("sintax_cutoff", 0.8))
        db_sources = set(
            (run_config_data or {}).get("db_sources")
            or (run_config_data or {}).get("vsearch_params", {}).get("db_sources")
            or []
        )
       
        silva_config = SILVAConfig(db_path=db_paths.silva) if "silva" in db_sources and db_paths.silva and db_paths.silva.exists() else None
        pr2_config = PR2Config(db_path=db_paths.pr2, tax_path=db_paths.pr2_tax) if "pr2" in db_sources and db_paths.pr2 and db_paths.pr2.exists() else None
        # BOLDigger3 works against the remote BOLD API, so it must NOT be gated on
        # a local directory existing — only on the user enabling the BOLD source.
        # ``run_boldigger`` surfaces a clear message if credentials/tool are absent.
        bold_config = (
            BOLDConfig(
                db_path=db_paths.bold_dir or Path("databases/bold"),
                db_source=int(run_config.bold_params.get("db", 1)),
                mode=int(run_config.bold_params.get("mode", 2)),
            )
            if "bold" in db_sources
            else None
        )
        ncbi_config = (
            NCBIConfig(
                local_db_path=db_paths.blast,
                use_remote=not (db_paths.blast and db_paths.blast.exists()),
                evalue=float(run_config.vsearch_params.get("ncbi_evalue", 1e-10)),
                max_target_seqs=int(run_config.vsearch_params.get("ncbi_maxhits", 5)),
            )
            if "ncbi" in db_sources
            else None
        )

        log_lines = [{"ts": utc_now_iso(), "level": "RUN", "msg": "Starting taxonomy pipeline"}]
        taxonomy_result = run_taxonomy_pipeline(
			fasta_path,
			run_config.output_dir / "databases",
			silva_config=silva_config,
			pr2_config=pr2_config,
			ncbi_config=ncbi_config,
			bold_config=bold_config,
            custom_db_path=custom_db_path,
            custom_db_kind=custom_db_kind,
            abundance_matrix_path=Path(abundance_matrix_path) if abundance_matrix_path else None,
            sintax_cutoff=sintax_cutoff,
            consensus_method=str(run_config.vsearch_params.get("db_consensus") or "lca"),
            threads=int(run_config.vsearch_params.get("threads", 8)),
            # Cut SILVA/PR2 to the region this run amplified: it is what lets them run
            # locally at all on a small machine, and it matches reference to query.
            marker=resolve_marker(sample.marker for sample in run_config.samples),
            primers=resolve_primers(run_config),
		)
        log_lines.append({"ts": utc_now_iso(), "level": "OK" if taxonomy_result.success else "ERR", "msg": taxonomy_result.error_message or "Taxonomy pipeline completed"})
        try:
            from databases.database_metrics import build_database_metrics

            metrics_files = {
                key: str(path)
                for key, path in build_database_metrics(
                    run_config.output_dir / "databases",
                    taxonomy_result.files.get("taxonomy_consensus_csv"),
                    abundance_matrix_path,
                    control_samples=[sample.internal_id for sample in run_config.samples if sample.is_control],
                    groups={sample.internal_id: sample.group for sample in run_config.samples if sample.group},
                ).items()
            }
        except Exception as exc:
            metrics_files = {}
            log_lines.append({"ts": utc_now_iso(), "level": "ERR", "msg": f"Falha ao construir métricas de taxonomia: {exc}"})

        record = {
            "signature": signature,
            "status": "success" if taxonomy_result.success else "failure",
            "updated_at": utc_now_iso(),
            "outputs": {"taxonomy": serialize_step_output(taxonomy_result)},
            "metrics_files": metrics_files,
            "log_lines": log_lines,
            "error_message": taxonomy_result.error_message,
        }
        return run_db_classification(project_data=project_data, samples_data=samples_data, step_outputs={"analysis": analysis_record, "databases": record}), record
