"""Results-page summary callbacks: run header, read funnel, stage cards, exports.

Everything here reads the run from disk (core.run_registry), not from the session
stores, so Results is meaningful after a restart or on a fresh tab.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from ui.components.cards import empty_state, kpi, status_badge
from ui.components.table import data_table

_GRAPH_CONFIG = {"displaylogo": False, "responsive": True}


def _state(run_config_data, project_data):
    from core.run_registry import load_run_state, run_dir_for

    run_dir = run_dir_for(run_config_data, project_data)
    return run_dir, load_run_state(run_dir)


def _read_csv(path):
    import pandas as pd

    try:
        return pd.read_csv(path)
    except Exception:
        return None


def _metric_paths(state: dict, stage: str) -> dict:
    from core.run_registry import stage_record

    files = stage_record(state, stage).get("metrics_files") or {}
    return {key: value for key, value in files.items() if isinstance(value, str) and Path(value).exists()}


def register_callbacks(app: Any) -> None:

    # ── Bloco A — run header ─────────────────────────────────────────────
    @app.callback(
        Output("results-run-header", "children"),
        Input("store-run-config", "data"),
        Input("store-databases-outputs", "data"),
        State("store-project", "data"),
    )
    def _render_header(run_config_data, _db, project_data):
        from core.run_registry import overall_status, run_stages

        run_dir, state = _state(run_config_data, project_data)
        if not state:
            return empty_state("Nenhum run encontrado.", hint="Corre o pipeline em «Set up your run».")

        overall = overall_status(state)
        stages = run_stages(state)
        run_id = state.get("run_id") or (run_dir.name if run_dir else "—")
        project = run_dir.parent.name if run_dir else "—"
        db_sources = (run_config_data or {}).get("db_sources") or []

        chips = [
            html.Div(f"Projeto: {project}", className="status-pill"),
            html.Div(f"Run: {run_id}", className="status-pill"),
            html.Div(f"Bases: {', '.join(db_sources).upper() or 'nenhuma'}", className="status-pill"),
        ]
        summary = " · ".join(f"{s['label']} {s['status']}" for s in stages)
        return html.Div(
            [
                html.Div(
                    [
                        html.Div("Estado do run", className="card-title", style={"marginBottom": "0"}),
                        status_badge(overall),
                    ],
                    style={"display": "flex", "alignItems": "center", "gap": "12px"},
                ),
                html.Div(chips, className="chip-list", style={"marginTop": "8px"}),
                html.Div(summary, className="mini-note"),
            ]
        )

    # ── Bloco B — read funnel ────────────────────────────────────────────
    @app.callback(
        Output("results-funnel", "children"),
        Input("store-run-config", "data"),
        Input("store-databases-outputs", "data"),
        State("store-project", "data"),
    )
    def _render_funnel(run_config_data, _db, project_data):
        from ui.plots import read_funnel_figure

        run_dir, state = _state(run_config_data, project_data)
        if not state:
            return empty_state("Sem run.")

        steps = _funnel_steps(run_dir, state)
        if not steps:
            return empty_state("Sem dados de tracking suficientes para construir o funil.")
        return dcc.Graph(figure=read_funnel_figure(steps), config=_GRAPH_CONFIG)

    # ── Bloco C — per-stage KPI cards ────────────────────────────────────
    @app.callback(
        Output("results-stage-cards", "children"),
        Input("store-run-config", "data"),
        Input("store-databases-outputs", "data"),
        State("store-project", "data"),
    )
    def _render_stage_cards(run_config_data, _db, project_data):
        from core.run_registry import run_stages, stage_record

        run_dir, state = _state(run_config_data, project_data)
        if not state:
            return empty_state("Sem run.")

        by_key = {stage["key"]: stage for stage in run_stages(state)}
        cards = []

        for key, href, builder in (
            ("qc", "/qc", _qc_card),
            ("preprocess", "/preprocess", _preprocess_card),
            ("analysis", "/analysis", _analysis_card),
            ("databases", "/databases", _databases_card),
            ("phylo", "/phylo", _phylo_card),
        ):
            stage = by_key.get(key, {})
            status = stage.get("status", "pending")
            if status in ("skipped", "pending") or not stage.get("has_outputs"):
                # A stage that did not run says so. Never a blank card.
                cards.append(
                    html.Div(
                        [
                            html.Div("—", className="stat-value", style={"color": "#8b98aa"}),
                            html.Div(stage.get("label", key), className="stat-label"),
                            html.Div(status_badge(status), style={"marginTop": "4px"}),
                        ],
                        className="stat-card",
                    )
                )
                continue
            cards.append(builder(run_dir, state, href, status))

        return html.Div(cards, className="grid-5" if len(cards) == 5 else "grid-3")

    # ── Bloco D — composition ────────────────────────────────────────────
    @app.callback(
        Output("results-composition", "children"),
        Input("results-composition-rank", "value"),
        Input("results-composition-exclude", "value"),
        Input("store-databases-outputs", "data"),
        State("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_composition(rank, exclude, _db, run_config_data, project_data):
        from ui.plots import stacked_composition_figure

        run_dir, state = _state(run_config_data, project_data)
        if not run_dir:
            return empty_state("Sem run.")

        frame = _composition_table(run_dir, rank or "Consensus_Phylum", bool(exclude))
        if frame is None or frame.empty:
            return empty_state(
                "Sem taxonomia. Corre a etapa Databases com pelo menos uma base de referência.",
            )
        return dcc.Graph(figure=stacked_composition_figure(frame, rank or "Consensus_Phylum"), config=_GRAPH_CONFIG)

    # ── Bloco G — exports ────────────────────────────────────────────────
    @app.callback(
        Output("results-exports", "children"),
        Input("store-run-config", "data"),
        Input("store-databases-outputs", "data"),
        State("store-project", "data"),
    )
    def _render_exports(run_config_data, _db, project_data):
        run_dir, state = _state(run_config_data, project_data)
        if not run_dir:
            return empty_state("Sem run.")

        from core.config import get_runs_dir

        runs_root = Path(get_runs_dir())
        rows = []
        for label, relative, reason in _EXPORTS:
            path = run_dir / relative
            if path.exists():
                size = path.stat().st_size
                try:
                    href = f"/runs/{path.relative_to(runs_root).as_posix()}"
                except ValueError:
                    href = None
                rows.append(
                    {
                        "artefacto": label,
                        "ficheiro": path.name,
                        "tamanho": _human_size(size),
                        "estado": "disponível",
                        "_href": href,
                    }
                )
            else:
                # Missing artifacts are listed greyed-out WITH the reason, rather than
                # omitted — otherwise the user cannot tell "not produced" from "bug".
                rows.append({"artefacto": label, "ficheiro": "—", "tamanho": "—", "estado": reason, "_href": None})

        links = [
            html.Div(html.A(f"↓ {row['artefacto']} — {row['ficheiro']} ({row['tamanho']})", href=row["_href"], target="_blank"), className="download-row")
            for row in rows
            if row["_href"]
        ]
        table = data_table(
            data=[{k: v for k, v in row.items() if not k.startswith("_")} for row in rows],
            columns=[{"name": n.capitalize(), "id": n} for n in ("artefacto", "ficheiro", "tamanho", "estado")],
            conditional=[
                {"if": {"filter_query": '{estado} != "disponível"'}, "color": "#8b98aa", "fontStyle": "italic"},
            ],
        )
        return html.Div([table, html.Div(links, style={"marginTop": "10px"})])

    @app.callback(
        Output("results-download-all", "data"),
        Input("btn-download-all", "n_clicks"),
        State("store-run-config", "data"),
        State("store-project", "data"),
        prevent_initial_call=True,
    )
    def _download_all(n_clicks, run_config_data, project_data):
        if not n_clicks:
            raise PreventUpdate
        run_dir, _ = _state(run_config_data, project_data)
        if not run_dir:
            raise PreventUpdate
        return _zip_download(run_dir)

    # ── Bloco H — reproducibility ────────────────────────────────────────
    @app.callback(
        Output("results-repro", "children"),
        Input("store-run-config", "data"),
        State("store-project", "data"),
    )
    def _render_repro(run_config_data, project_data):
        run_dir, state = _state(run_config_data, project_data)
        if not state:
            return empty_state("Sem run.")
        text = _reproducibility_block(run_dir, state, run_config_data)
        return html.Div(
            [
                html.Div("Copia este bloco para a secção de métodos.", className="mini-note"),
                dcc.Textarea(value=text, readOnly=True, style={"width": "100%", "height": "280px", "fontFamily": "monospace", "fontSize": "12px", "background": "#0e141c", "color": "#c2ccd8", "border": "1px solid #202b3c", "borderRadius": "8px", "padding": "10px"}),
            ]
        )


# ── helpers ──────────────────────────────────────────────────────────────

def _zip_download(run_dir: Path):
    import io

    from core.export import build_run_zip
    from core.run_registry import load_run_state

    state = load_run_state(run_dir)
    stage_outputs = {
        "qc": state.get("qc_record"),
        "preprocess": state.get("preprocess_record"),
        "analysis": state.get("analysis_record"),
        "databases": state.get("db_record"),
        "phylo": state.get("phylo_record"),
    }
    buf: io.BytesIO = build_run_zip(run_dir.name, stage_outputs)
    return dcc.send_bytes(buf.getvalue(), f"{run_dir.name}.zip")


def _human_size(size: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} GB"


# label, path relative to the run dir, reason it might be missing
#
# The taxonomy-derived rows below (relative composition, BIOM, phyloseq) need a
# consensus, not merely a Databases stage that ran: with every reference skipped for
# lack of RAM there is no taxonomy to collapse, and they are legitimately absent. Their
# reason says so, because "Databases não correu" would be a lie the user would go and
# try to fix.
_NO_TAXONOMY = "sem taxonomia — a etapa Databases não correu ou nenhuma fonte deu resultados"

_EXPORTS = [
    ("Tabela de OTUs", "analysis/otu_matrix.csv", "a etapa Analysis não correu"),
    ("Representativas (FASTA)", "analysis/pooled/otus.fasta", "a etapa Analysis não correu"),
    ("ASVs (FASTA)", "analysis/dada2/asv_sequences.fasta", "o DADA2 não correu"),
    ("Tabela de ASVs", "analysis/dada2/asv_abundance.csv", "o DADA2 não correu"),
    ("Taxonomia (consenso)", "databases/otu/taxonomy_consensus.csv", "a etapa Databases não correu"),
    ("Todos os hits", "databases/otu/taxonomy_hits.csv", "a etapa Databases não correu"),
    ("Diversidade alfa", "analysis/metrics/otu/alpha_diversity.csv", "a etapa Analysis não correu"),
    ("Curvas de rarefação", "analysis/metrics/otu/rarefaction.csv", "a etapa Analysis não correu"),
    ("Composição relativa (género, %)", "databases/otu/taxa_tables/taxa_table_genus_relative.csv", _NO_TAXONOMY),
    ("Composição relativa (espécie, %)", "databases/otu/taxa_tables/taxa_table_species_relative.csv", _NO_TAXONOMY),
    ("Tabela BIOM (QIIME/phyloseq)", "databases/otu/table.biom", _NO_TAXONOMY),
    # The link gives the script alone; its three tables are siblings on disk and travel
    # in the run .zip. The label says "script" rather than "script + tabelas" because
    # the download is the script — promising the tables in a single click would be a lie.
    ("phyloseq (script R, tabelas no .zip)", "databases/otu/phyloseq/load_phyloseq.R", _NO_TAXONOMY),
    # The phylo artifacts sit under their runner's subdir (MAFFT/RAxML), not flat in the
    # set dir — the previous ``phylo/otu/tree.nwk`` never existed, so a completed tree
    # always read as "não correu". These paths match the default HYBRID mode the golden
    # freezes; the support tree is the Newick the viewer surfaces.
    ("Alinhamento", "phylo/otu/mafft/aligned.fasta", "a etapa Phylo não correu"),
    ("Árvore (Newick, com suporte)", "phylo/otu/raxml/raxml.raxml.support", "a etapa Phylo não correu"),
]


def _funnel_steps(run_dir: Path, state: dict) -> list[dict]:
    """The read funnel, assembled from the tables each stage already persisted."""
    from core.run_registry import stage_record

    steps: list[dict] = []

    tracking = _metric_paths(state, "preprocess").get("read_tracking_csv")
    frame = _read_csv(tracking) if tracking else None
    if frame is not None and not frame.empty:
        for column, label in (
            ("input_reads", "Input"),
            ("trimmed_reads", "Após trimming"),
            ("assembled_reads", "Merged (PEAR)"),
        ):
            if column in frame.columns:
                total = float(frame[column].fillna(0).sum())
                if total > 0:
                    steps.append({"label": label, "reads": total})

    # Clustering and classification totals come from the analysis/databases artifacts.
    from analysis.results_summary import clustered_read_total, species_classified_total

    clustered = clustered_read_total(run_dir)
    if clustered:
        steps.append({"label": "Clustered (OTU/ASV)", "reads": float(clustered)})
    classified = species_classified_total(run_dir)
    if classified:
        steps.append({"label": "Classificadas", "reads": float(classified)})

    return steps if len(steps) >= 2 else []


def _composition_table(run_dir: Path, rank: str, exclude_host: bool):
    """Relative abundance per sample at ``rank``, as a taxon x sample frame."""
    import pandas as pd

    from analysis.analysis_metrics import load_matrix

    consensus = run_dir / "databases" / "otu" / "taxonomy_consensus.csv"
    if not consensus.exists():
        consensus = run_dir / "databases" / "taxonomy_consensus.csv"
    matrix_path = run_dir / "analysis" / "otu_matrix.csv"
    if not consensus.exists() or not matrix_path.exists():
        return None

    taxonomy = _read_csv(consensus)
    if taxonomy is None or rank not in taxonomy.columns:
        return None

    features, samples, counts = load_matrix(matrix_path)
    if not features or not samples:
        return None

    id_col = "feature_id" if "feature_id" in taxonomy.columns else "ASV_ID"
    if id_col not in taxonomy.columns:
        return None
    names = dict(zip(taxonomy[id_col].astype(str), taxonomy[rank]))
    host = _HOST_MARKERS

    rows = []
    for index, feature in enumerate(features):
        # A rank column with no assignments comes back as float NaN, not "", so coerce
        # before touching it — .strip() on a float is what took this callback down.
        taxon = str(names.get(feature, "") or "").strip()
        if not taxon or taxon.lower() in ("nan", "none"):
            taxon = "Não classificado"
        if exclude_host:
            lineage = taxon.lower()
            if taxon == "Não classificado" or any(marker in lineage for marker in host):
                continue
        rows.append({"taxon": taxon, **{sample: counts[s, index] for s, sample in enumerate(samples)}})

    if not rows:
        return None
    frame = pd.DataFrame(rows).groupby("taxon", as_index=False).sum()
    for sample in samples:
        total = frame[sample].sum()
        frame[sample] = frame[sample] / total * 100 if total else 0
    return frame


_HOST_MARKERS = ("chordata", "craniata", "actinopteri", "vertebrata", "mammalia", "micromesistius")


# Each tool announces its version differently, and several exit non-zero or print a
# citation when asked. A per-tool flag + pattern is the only way to get a version rather
# than "USAGE" or a bibliography entry into someone's methods section.
#
# The executable is resolved from core.config, NOT from $PATH: the pipeline runs the
# binary that BLASTN_PATH & friends point at, and reporting whatever `which` happens to
# find would put a version in a paper that never touched the data. On this machine those
# differ — $PATH has the broken blastn 2.9.0, while the pipeline uses 2.17.
_TOOL_VERSIONS = [
    ("vsearch", "VSEARCH_PATH", "--version", r"vsearch\s+v?([\d.]+)"),
    ("blastn", "BLASTN_PATH", "-version", r"blastn:\s*([\d.+]+)"),
    ("fastqc", "FASTQC_PATH", "--version", r"FastQC\s+v?([\d.]+)"),
    ("mafft", "MAFFT_PATH", "--version", r"v([\d.]+)"),
    ("FastTree", "FASTTREE_PATH", "", r"(?:FastTree|Detected)\s+(?:version\s+)?([\d.]+)"),
    ("RAxML-NG", "RAXML_PATH", "--version", r"RAxML-NG v\.\s*([\d.]+)"),
    ("Rscript", "RSCRIPT_PATH", "--version", r"version\s+([\d.]+)"),
    ("PEAR", "PEAR_PATH", "--version", r"PEAR\s+v?([\d.]+)"),
]


def _reproducibility_block(run_dir: Path, state: dict, run_config_data: Optional[dict]) -> str:
    """A methods-section paragraph the user can paste into a paper."""
    import re
    import shutil
    import subprocess

    config = run_config_data or {}
    lines = [
        "PARAFISH — bloco de reprodutibilidade",
        "=" * 46,
        f"Run ID       : {state.get('run_id', '—')}",
        f"Directório   : {run_dir}",
        f"Método       : {config.get('analysis_method', '—')}",
        f"Bases        : {', '.join(config.get('db_sources') or []) or '—'}",
        "",
        "Ferramentas",
        "-" * 46,
    ]

    import core.config as tool_config

    for tool, setting, flag, pattern in _TOOL_VERSIONS:
        path = getattr(tool_config, setting, "") or ""
        resolved = shutil.which(path) or (path if Path(path).is_file() else None)
        if not resolved:
            lines.append(f"  {tool:10s} não instalado")
            continue
        try:
            # FastTree has no version flag: it prints it in its banner on a bare call.
            argv = [resolved, flag] if flag else [resolved]
            out = subprocess.run(argv, capture_output=True, text=True, timeout=15, input="")
            stream = f"{out.stdout or ''}\n{out.stderr or ''}"
            match = re.search(pattern, stream)
            version = match.group(1) if match else "instalado (versão não detetada)"
        except Exception:
            version = "instalado (versão não detetada)"
        lines.append(f"  {tool:10s} {version:12s} {resolved}")

    lines += ["", "Parâmetros", "-" * 46]
    for section in ("trimmomatic_params", "pear_params", "vsearch_params", "dada2_params", "phylo_params"):
        params = config.get(section)
        if isinstance(params, dict) and params:
            lines.append(f"  [{section}]")
            for key, value in params.items():
                if key in ("db_sources", "custom_db_path"):
                    continue
                lines.append(f"    {key} = {value}")

    lines += ["", "Ordenação", "-" * 46]
    from analysis.ordination import MIN_PREVALENCE, ZERO_REPLACEMENT

    lines.append(f"  Método por omissão : PCA de Aitchison (CLR)")
    lines.append(f"  Substituição de zeros : pseudocount = {ZERO_REPLACEMENT} (sobre contagens)")
    lines.append(f"  Prevalência mínima : {MIN_PREVALENCE} amostras")

    depth = _metric_paths(state, "analysis").get("rarefaction_depth_csv") or _metric_paths(state, "analysis").get("otu.rarefaction_depth_csv")
    frame = _read_csv(depth) if depth else None
    if frame is not None and not frame.empty and "depth" in frame.columns:
        lines.append(f"  Profundidade de rarefação : {int(frame['depth'].iloc[0]):,} reads")

    return "\n".join(lines)


# ── per-stage KPI card builders ──────────────────────────────────────────

def _qc_card(run_dir, state, href, status):
    paths = _metric_paths(state, "qc")
    overview = _read_csv(paths.get("qc_overview_csv")) if paths.get("qc_overview_csv") else None
    value, hint = "—", None
    if overview is not None and not overview.empty:
        modules = _read_csv(paths.get("qc_module_status_csv")) if paths.get("qc_module_status_csv") else None
        if modules is not None and "status" in modules.columns:
            bad = int(modules["status"].astype(str).str.lower().isin(["warn", "fail"]).sum())
            total = len(modules)
            value = f"{(1 - bad / total) * 100:.0f}%" if total else "—"
            hint = f"{bad} módulos em warn/fail"
        else:
            value = f"{len(overview)} ficheiros"
    return kpi(value, "QC — pass rate", hint=hint, href=href)


def _preprocess_card(run_dir, state, href, status):
    paths = _metric_paths(state, "preprocess")
    frame = _read_csv(paths.get("read_tracking_csv")) if paths.get("read_tracking_csv") else None
    value, hint = "—", None
    if frame is not None and not frame.empty and "pct_surviving" in frame.columns:
        value = f"{frame['pct_surviving'].mean():.0f}%"
        if "pct_assembled" in frame.columns:
            hint = f"assembly médio {frame['pct_assembled'].mean():.0f}%"
    return kpi(value, "Preprocess — reads retidas", hint=hint, href=href)


def _analysis_card(run_dir, state, href, status):
    from core.run_registry import stage_record

    outputs = stage_record(state, "analysis").get("outputs") or {}
    sets = outputs.get("feature_sets") or []
    counts = []
    for entry in sets:
        fasta = Path(entry.get("fasta", ""))
        if fasta.exists():
            n = sum(1 for line in fasta.read_text(errors="replace").splitlines() if line.startswith(">"))
            counts.append(f"{n} {entry['key'].upper()}s")
    value = " · ".join(counts) if counts else "—"

    paths = _metric_paths(state, "analysis")
    chimera = _read_csv(paths.get("otu.chimera_stats_csv") or paths.get("chimera_stats_csv")) if paths else None
    hint = None
    if chimera is not None and not chimera.empty:
        for column in ("pct_chimeric_reads", "chimeric_read_pct", "pct_reads"):
            if column in chimera.columns:
                hint = f"{float(chimera[column].iloc[0]):.1f}% reads quiméricas"
                break
    return kpi(value, "Analysis — features", hint=hint, href=href)


def _databases_card(run_dir, state, href, status):
    consensus = run_dir / "databases" / "otu" / "taxonomy_consensus.csv"
    if not consensus.exists():
        consensus = run_dir / "databases" / "taxonomy_consensus.csv"
    frame = _read_csv(consensus) if consensus.exists() else None
    value, hint = "—", None
    if frame is not None and not frame.empty and "Consensus_Genus" in frame.columns:
        genus = frame["Consensus_Genus"].astype(str).str.strip()
        named = (genus != "") & (genus.str.lower() != "nan")
        value = f"{int(named.sum()) / len(frame) * 100:.0f}%"
        conflicts = run_dir / "databases" / "otu" / "conflict_report.csv"
        if conflicts.exists():
            report = _read_csv(conflicts)
            if report is not None:
                hint = f"{len(report)} conflitos entre fontes"
    return kpi(value, "Databases — classificadas a género", hint=hint, href=href)


def _phylo_card(run_dir, state, href, status):
    paths = _metric_paths(state, "phylo")
    value, hint = "—", None
    from ui.feature_sets import phylo_tree_path

    support = _read_csv(paths.get("otu.support_csv") or paths.get("support_csv")) if paths else None
    tree = phylo_tree_path(run_dir, "otu")
    if tree:
        text = Path(tree).read_text(errors="replace")
        value = f"{text.count(',') + 1} tips"
    if support is not None and not support.empty and "support" in support.columns:
        hint = f"suporte mediano {support['support'].median():.0f}"
    return kpi(value, "Phylo — sequências na árvore", hint=hint, href=href)
