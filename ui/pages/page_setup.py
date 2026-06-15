"""Unified pipeline setup page."""

from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html


def _param_label(text: str, tooltip: str):
    return html.Div(
        [
            html.Label(text, className="form-label"),
            html.Span(
                "ⓘ",
                title=tooltip,
                style={
                    "cursor": "help",
                    "marginLeft": "6px",
                    "color": "var(--text-2)",
                    "fontSize": "12px",
                },
            ),
        ],
        style={"display": "flex", "alignItems": "center"},
    )


def _workflow_mode_card():
    return html.Div(
        [
            html.Div("Workflow mode", className="card-title"),
            dcc.RadioItems(
                id="setup-workflow-mode",
                options=[
                    {"label": "Complete", "value": "complete"},
                    {"label": "Special", "value": "special"},
                    {"label": "One", "value": "one"},
                ],
                value="complete",
                className="radio-group",
                labelStyle={"display": "inline-block", "marginRight": "12px"},
            ),
        ],
        className="card card-body",
        style={"marginBottom": "12px"},
    )


def _qc_card():
    return html.Div(
        [
            html.Div(
                [
                    html.Span("QC (FastQC)", className="card-title"),
                    dbc.Switch(
                        id="setup-enable-qc",
                        label="ON",
                        value=True,
                    ),
                ],
                className="card-header",
            ),
            html.Div("FastQC has no manual parameters.", className="empty-state-sub card-body"),
        ],
        className="card",
        style={"marginBottom": "12px"},
    )


def _preprocess_card():
    return html.Div(
        [
            html.Div(
                [
                    html.Span("Preprocess", className="card-title"),
                    dbc.Switch(
                        id="setup-enable-preprocess",
                        label="ON",
                        value=True,
                    ),
                ],
                className="card-header",
            ),
            html.Div(
                [
                    _param_label("LEADING", "Removes low-quality bases from the 5' end. Higher = more aggressive. Default: 3."),
                    dcc.Input(id="setup-trim-leading", type="number", value=3, min=0, max=30, className="dash-input"),
                    _param_label("TRAILING", "Removes low-quality bases from the 3' end. Higher = more aggressive. Default: 3."),
                    dcc.Input(id="setup-trim-trailing", type="number", value=3, min=0, max=30, className="dash-input"),
                    _param_label("SLIDINGWINDOW", "Window:quality — scans a window of N bases and cuts when average quality drops below threshold. Default: 4:15."),
                    dcc.Input(id="setup-trim-slidingwindow", type="text", value="4:15", className="dash-input"),
                    _param_label("MINLEN", "Minimum read length after trimming. Reads shorter than this are discarded. Default: 50."),
                    dcc.Input(id="setup-trim-minlen", type="number", value=50, min=20, max=500, className="dash-input"),
                    html.Hr(),
                    html.Div("PEAR (PE only)", className="param-section-title"),
                    _param_label("Min overlap", "(PEAR) Minimum overlap between forward and reverse reads for merging. Default: 10."),
                    dcc.Input(id="setup-pear-min-overlap", type="number", value=10, min=1, max=100, className="dash-input"),
                    _param_label("Min length", "(PEAR) Minimum length of the merged read. Default: 200."),
                    dcc.Input(id="setup-pear-min-length", type="number", value=200, min=20, max=1000, className="dash-input"),
                    _param_label("P-value", "(PEAR) Statistical threshold for assembly acceptance. Lower = stricter. Default: 0.01."),
                    dcc.Input(id="setup-pear-pvalue", type="number", value=0.01, min=0.0, max=1.0, step=0.001, className="dash-input"),
                ],
                className="card-body",
            ),
        ],
        className="card",
        style={"marginBottom": "12px"},
    )


def _analysis_card():
    return html.Div(
        [
            html.Div(
                [
                    html.Span("Analysis", className="card-title"),
                    dbc.Switch(
                        id="setup-enable-analysis",
                        label="ON",
                        value=True,
                    ),
                ],
                className="card-header",
            ),
            html.Div(
                [
                    html.Label("Method", className="form-label"),
                    dcc.RadioItems(
                        id="setup-analysis-method",
                        options=[
                            {"label": "VSEARCH (OTUs)", "value": "vsearch"},
                            {"label": "DADA2 (ASVs)", "value": "dada2"},
                            {"label": "Both", "value": "both"},
                        ],
                        value="both",
                        className="radio-group",
                        labelStyle={"display": "block", "marginBottom": "6px"},
                    ),
                    html.Div(
                        id="vsearch-params-panel",
                        children=[
                            html.Div("VSEARCH", className="param-section-title"),
                            _param_label("Identity threshold (--id)", "(VSEARCH) Minimum sequence similarity to cluster reads into an OTU. 0.97 = species-level OTU."),
                            dcc.Slider(
                                id="setup-vsearch-id",
                                min=0.80,
                                max=1.00,
                                step=0.01,
                                value=0.97,
                                marks={0.80: "0.80", 0.97: "0.97 OTU", 0.99: "0.99 ASV", 1.00: "1.00"},
                                tooltip={"placement": "bottom", "always_visible": True},
                            ),
                            _param_label("Max expected errors", "(VSEARCH) Maximum number of expected errors per read. Higher = more permissive. Default: 1.0."),
                            dcc.Input(id="setup-vsearch-maxee", type="number", value=1.0, min=0.1, max=5.0, step=0.1, className="dash-input"),
                            _param_label("Min length (bp)", "(VSEARCH) Minimum amplicon length to retain after quality filtering. Default: 200."),
                            dcc.Input(id="setup-vsearch-minlen", type="number", value=200, min=50, max=600, step=10, className="dash-input"),
                            _param_label("SINTAX cutoff", "(VSEARCH) Minimum bootstrap confidence to assign taxonomy. 0.8 = 80% confidence required."),
                            dcc.Input(id="setup-vsearch-sintax", type="number", value=0.8, min=0.5, max=1.0, step=0.05, className="dash-input"),
                        ],
                        style={"display": "block", "marginTop": "12px"},
                    ),
                    html.Div(
                        id="dada2-params-panel",
                        children=[
                            html.Div("DADA2", className="param-section-title"),
                            _param_label("TruncLen Forward", "(DADA2) Truncate forward reads to this length. Set 0 to disable. Disabled for ITS."),
                            dcc.Input(id="setup-dada2-truncF", type="number", value=0, min=0, max=300, className="dash-input"),
                            _param_label("TruncLen Reverse", "(DADA2) Truncate reverse reads to this length. Set 0 to disable."),
                            dcc.Input(id="setup-dada2-truncR", type="number", value=0, min=0, max=300, className="dash-input"),
                            _param_label("MaxEE Forward", "(DADA2) Maximum expected errors in forward read. Reads above this are discarded. Default: 2."),
                            dcc.Input(id="setup-dada2-maxeeF", type="number", value=2, min=1, max=10, className="dash-input"),
                            _param_label("MaxEE Reverse", "(DADA2) Maximum expected errors in reverse read. Default: 2."),
                            dcc.Input(id="setup-dada2-maxeeR", type="number", value=2, min=1, max=10, className="dash-input"),
                            html.Div(
                                id="dada2-its-warning",
                                children="ITS marker detected: truncation disabled.",
                                className="badge badge-amber",
                                style={"display": "none", "marginTop": "8px"},
                            ),
                        ],
                        style={"display": "block", "marginTop": "12px"},
                    ),
                ],
                className="card-body",
            ),
        ],
        className="card",
        style={"marginBottom": "12px"},
    )


def _databases_card():
    return html.Div(
        [
            html.Div(
                [
                    html.Span("Databases", className="card-title"),
                    dbc.Switch(
                        id="setup-enable-databases",
                        label="ON",
                        value=True,
                    ),
                ],
                className="card-header",
            ),
            html.Div(
                [
                    _param_label("Sources", "Database sources to include in the taxonomy step. Choose based on your markers."),
                    html.Div(id="setup-db-recommendation", className="badge badge-blue", style={"marginBottom": "8px"}),
                    dcc.Checklist(
                        id="setup-db-sources",
                        options=[
                            {"label": "SILVA", "value": "silva"},
                            {"label": "PR2", "value": "pr2"},
                            {"label": "BOLD", "value": "bold"},
                            {"label": "NCBI", "value": "ncbi"},
                        ],
                        value=["silva", "pr2", "bold"],
                        className="checkbox-group",
                        labelStyle={"display": "inline-block", "marginRight": "10px"},
                    ),
                    html.Div(
                        id="db-params-container",
                        children=[
                            html.Div(
                                [
                                    html.Label("SINTAX confidence cutoff (SILVA/PR2)"),
                                    dcc.Slider(
                                        id="setup-sintax-cutoff",
                                        min=0.5,
                                        max=1.0,
                                        step=0.05,
                                        value=0.8,
                                        marks={0.5: "50%", 0.8: "80%", 1.0: "100%"},
                                        tooltip={"placement": "bottom", "always_visible": True},
                                    ),
                                    html.Label("NCBI e-value"),
                                    dcc.Input(id="setup-ncbi-evalue", type="text", value="1e-10", style={"width": "120px"}),
                                    html.Label("NCBI max hits"),
                                    dcc.Input(id="setup-ncbi-maxhits", type="number", value=5, min=1, max=50, style={"width": "80px"}),
                                    html.Label("BOLD database number (1-8)"),
                                    dcc.Dropdown(
                                        id="setup-bold-db",
                                        options=[
                                            {"label": "1 — Animal Library (Public)", "value": 1},
                                            {"label": "2 — Animal Species-level (Public+Private)", "value": 2},
                                            {"label": "5 — Plant Library (Public)", "value": 5},
                                            {"label": "6 — Fungi Library (Public)", "value": 6},
                                        ],
                                        value=1,
                                        clearable=False,
                                    ),
                                    html.Label("BOLD mode"),
                                    dcc.Dropdown(
                                        id="setup-bold-mode",
                                        options=[
                                            {"label": "1 — Rapid Species Search", "value": 1},
                                            {"label": "2 — Genus and Species Search", "value": 2},
                                            {"label": "3 — Exhaustive Search", "value": 3},
                                        ],
                                        value=1,
                                        clearable=False,
                                    ),
                                    html.Label("Consensus method"),
                                    dcc.Dropdown(
                                        id="setup-db-consensus",
                                        options=[
                                            {"label": "LCA (Lowest Common Ancestor)", "value": "lca"},
                                            {"label": "Strict (all sources must agree)", "value": "strict"},
                                            {"label": "Majority vote", "value": "majority"},
                                        ],
                                        value="lca",
                                        clearable=False,
                                    ),
                                    _param_label("Species threshold (%)", "Minimum percent identity to assign at species level. Default: 97%."),
                                    dcc.Input(id="setup-species-threshold", type="number", value=97, min=50, max=100, step=1, className="dash-input"),
                                ],
                                style={"display": "grid", "gap": "8px"},
                            )
                        ],
                    ),
                ],
                className="card-body",
            ),
        ],
        className="card",
        style={"marginBottom": "12px"},
    )


def _phylo_card():
    return html.Div(
        [
            html.Div(
                [
                    html.Span("Phylogenetics", className="card-title"),
                    dbc.Switch(
                        id="setup-enable-phylo",
                        label="ON",
                        value=True,
                    ),
                ],
                className="card-header",
            ),
            html.Div(
                [
                    _param_label("Mode", "Choose the phylogenetic pipeline strategy: fast, robust, or hybrid."),
                    dcc.RadioItems(
                        id="setup-phylo-mode",
                        options=[
                            {"label": "Fast", "value": "fast"},
                            {"label": "Robust", "value": "robust"},
                            {"label": "Hybrid", "value": "hybrid"},
                        ],
                        value="hybrid",
                        className="radio-group",
                        labelStyle={"display": "inline-block", "marginRight": "10px"},
                    ),
                    _param_label("MAFFT", "Alignment strategy used before phylogenetic tree inference."),
                    dcc.RadioItems(
                        id="setup-mafft-mode",
                        options=[
                            {"label": "Auto", "value": "auto"},
                            {"label": "Accurate", "value": "accurate"},
                        ],
                        value="auto",
                        className="radio-group",
                        labelStyle={"display": "inline-block", "marginRight": "10px"},
                    ),
                    _param_label("RAxML Bootstraps", "Number of bootstrap replicates used to assess branch support. Higher values improve precision but increase runtime. Recommended: 10-50 for testing, 100-500 for final analyses."),
                    dcc.Input(id="setup-raxml-bootstraps", type="number", min=10, max=1000, step=10, value=10, className="dash-input"),
                    dcc.Checklist(
                        id="setup-adjust-direction",
                        options=[{"label": "Adjust direction", "value": "on"}],
                        value=["on"],
                        className="checkbox-group",
                    ),
                ],
                className="card-body",
            ),
        ],
        className="card",
        style={"marginBottom": "12px"},
    )


def _run_button():
    return html.Div(
        [
            html.Div(
                [
                    _param_label("Threads", "Number of CPU threads to use per tool. Match to your machine's available cores."),
                    dcc.Input(id="setup-threads", type="number", value=8, min=1, max=128, className="dash-input"),
                    _param_label("Workers", "Number of parallel sample workers. Higher = more RAM usage. Recommended: 2–4."),
                    dcc.Input(id="setup-workers", type="number", value=4, min=1, max=64, className="dash-input"),
                ],
                className="control-row",
            ),
            html.Button("Run Pipeline", id="btn-run-pipeline", className="btn-run"),
        ],
        className="card card-body",
    )


def layout():
    return html.Div(
        [
            html.Div(
                [
                    html.Div("Pipeline Configuration", className="page-title"),
                    html.Div("Configure parameters and run the full pipeline", className="page-sub"),
                ],
                className="page-header",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            _workflow_mode_card(),
                            _qc_card(),
                            _preprocess_card(),
                            _analysis_card(),
                            _databases_card(),
                            _phylo_card(),
                            _run_button(),
                        ],
                        style={"width": "400px", "flexShrink": "0"},
                    ),
                    html.Div(
                        [
                            html.Div(id="setup-samples-table"),
                            html.Div(id="setup-sample-count", className="badge badge-blue", style={"marginTop": "8px"}),
                        ],
                        style={"flex": "1", "minWidth": "0"},
                    ),
                ],
                style={"display": "flex", "gap": "20px", "alignItems": "flex-start"},
            ),
            html.Div(id="setup-run-status"),
        ],
        className="main-content",
    )
