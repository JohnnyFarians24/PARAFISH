"""Unified pipeline setup page."""

from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html

from core.config import MIN_UNIQUE_SIZE
from core.models import AnalysisBranch, RunConfig

# RunConfig is the single source of every default a Setup control stands for. A literal
# typed into a control is a second source, and it is the one that reaches a UI run: the
# sliding window (4:15), the RAxML bootstraps (100) and the BOLD mode (1) each disagreed
# with RunConfig (4:20, 10, 2), so a UI run and a headless run used different parameters.
_RUN_DEFAULTS = RunConfig(run_id="", samples=[], branch=AnalysisBranch.VSEARCH)
_DEFAULT_SLIDINGWINDOW = f"{_RUN_DEFAULTS.trimmomatic_params['window_size']}:{_RUN_DEFAULTS.trimmomatic_params['window_quality']}"


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


def _advanced(children, label: str = "Advanced settings (optional)"):
    """Collapse rarely-changed parameters behind a native disclosure.

    Keeps every input (and its id) in the DOM so callbacks still read them, while
    decluttering the default view — sensible defaults cover the common case.
    """
    return html.Details(
        [
            html.Summary(label, className="advanced-summary"),
            html.Div(children, className="advanced-body"),
        ],
        className="advanced-block",
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
                    html.Span("1 · Check read quality", className="card-title"),
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
                    html.Span("2 · Clean & merge reads", className="card-title"),
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
                    html.Div("Trims adapters/low-quality bases (Trimmomatic) and merges paired reads (PEAR). Defaults work for most Illumina runs.", className="empty-state-sub"),
                    _advanced(
                        [
                            html.Div("Trimmomatic", className="param-section-title"),
                            _param_label("LEADING", "Removes low-quality bases from the 5' end. Higher = more aggressive. Default: 3."),
                            dcc.Input(id="setup-trim-leading", type="number", value=3, min=0, max=30, className="dash-input"),
                            _param_label("TRAILING", "Removes low-quality bases from the 3' end. Higher = more aggressive. Default: 3."),
                            dcc.Input(id="setup-trim-trailing", type="number", value=3, min=0, max=30, className="dash-input"),
                            _param_label("SLIDINGWINDOW", f"Window:quality — scans a window of N bases and cuts when average quality drops below threshold. Default: {_DEFAULT_SLIDINGWINDOW}."),
                            dcc.Input(id="setup-trim-slidingwindow", type="text", value=_DEFAULT_SLIDINGWINDOW, className="dash-input"),
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
                        ]
                    ),
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
                    html.Span("3 · Find OTUs / ASVs", className="card-title"),
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
                            html.Div("VSEARCH (OTUs)", className="param-section-title"),
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
                            _advanced(
                                [
                                    _param_label("Max expected errors", "(VSEARCH) Maximum number of expected errors per read. Higher = more permissive. Default: 1.0."),
                                    dcc.Input(id="setup-vsearch-maxee", type="number", value=1.0, min=0.1, max=5.0, step=0.1, className="dash-input"),
                                    _param_label("Min length (bp)", "(VSEARCH) Minimum amplicon length to retain after quality filtering. Default: 200."),
                                    dcc.Input(id="setup-vsearch-minlen", type="number", value=200, min=50, max=600, step=10, className="dash-input"),
                                    _param_label("SINTAX cutoff", "(VSEARCH) Minimum bootstrap confidence to assign taxonomy. 0.8 = 80% confidence required."),
                                    dcc.Input(id="setup-vsearch-sintax", type="number", value=0.8, min=0.5, max=1.0, step=0.05, className="dash-input"),
                                    _param_label(
                                        "Min unique size",
                                        "(VSEARCH) Drop unique sequences seen fewer than this many times before clustering. "
                                        "The default is 1 — singletons are KEPT. A rare parasite is present as few reads, and "
                                        "a floor of 2 deletes exactly that profile before clustering ever sees it; measured on "
                                        "this data it discarded a quarter of the library. Raise it only to buy back clustering "
                                        "time on a large host-dominated study, and read unmapped_reads to see what it cost.",
                                    ),
                                    # Not a literal: this field is what a run actually uses, and a number typed
                                    # here silently outranked PARAFISH_MIN_UNIQUE_SIZE. The 2 that used to sit
                                    # here is why every UI run clustered with singletons dropped while the
                                    # project's documented, deliberate default was 1.
                                    dcc.Input(id="setup-vsearch-minsize", type="number", value=MIN_UNIQUE_SIZE, min=1, max=10, step=1, className="dash-input"),
                                ]
                            ),
                            _param_label(
                                "PCR primers",
                                "Trimmomatic removes sequencing adapters, NOT the PCR primers. Left in place, the primer's "
                                "degenerate positions are treated as biological signal and can split one organism across "
                                "several OTUs. Pick the pair you amplified with, or enter the sequences.",
                            ),
                            dcc.Dropdown(
                                id="setup-primer-preset",
                                options=[
                                    {"label": "Auto — the standard pair for this run's marker", "value": "auto"},
                                    {"label": "None — do not strip primers", "value": ""},
                                    {"label": "18S V4 — TAReuk454FWD1 / TAReukREV3", "value": "18S_V4_TAReuk"},
                                    {"label": "COI — Folmer (LCO1490 / HCO2198)", "value": "COI_Folmer"},
                                    {"label": "COI — Leray (mlCOIintF / jgHCO2198)", "value": "COI_Leray"},
                                    {"label": "COI — mlCOIintF / LoboR1 (lab)", "value": "COI_Leray_LoboR1"},
                                ],
                                # "auto", not "": the old default was the same value as
                                # "do not strip primers", so an untouched control meant
                                # primers were kept — and looked like a deliberate choice.
                                value="auto",
                                clearable=False,
                                className="dash-dropdown",
                            ),
                            _advanced(
                                [
                                    _param_label("Forward primer", "Overrides the preset. The first N bases are stripped from each merged read."),
                                    dcc.Input(id="setup-primer-fwd", type="text", placeholder="e.g. CCAGCASCYGCGGTAATTCC", debounce=True, className="dash-input"),
                                    _param_label("Reverse primer", "Overrides the preset. The last N bases are stripped from each merged read."),
                                    dcc.Input(id="setup-primer-rev", type="text", placeholder="e.g. ACTTTCGTTCTTGATYRA", debounce=True, className="dash-input"),
                                ]
                            ),
                        ],
                        style={"display": "block", "marginTop": "12px"},
                    ),
                    html.Div(
                        id="dada2-params-panel",
                        children=[
                            html.Div("DADA2 (ASVs)", className="param-section-title"),
                            _advanced(
                                [
                                    _param_label("TruncLen Forward", "(DADA2) Truncate forward reads to this length. Set 0 to disable. Disabled for ITS."),
                                    dcc.Input(id="setup-dada2-truncF", type="number", value=0, min=0, max=300, className="dash-input"),
                                    _param_label("TruncLen Reverse", "(DADA2) Truncate reverse reads to this length. Set 0 to disable."),
                                    dcc.Input(id="setup-dada2-truncR", type="number", value=0, min=0, max=300, className="dash-input"),
                                    _param_label("MaxEE Forward", "(DADA2) Maximum expected errors in forward read. Reads above this are discarded. Default: 2."),
                                    dcc.Input(id="setup-dada2-maxeeF", type="number", value=2, min=1, max=10, className="dash-input"),
                                    _param_label("MaxEE Reverse", "(DADA2) Maximum expected errors in reverse read. Default: 2."),
                                    dcc.Input(id="setup-dada2-maxeeR", type="number", value=2, min=1, max=10, className="dash-input"),
                                ]
                            ),
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
                    html.Span("4 · Identify species", className="card-title"),
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
                    html.Div(
                        "Match your sequences against reference databases to find out which species are present.",
                        className="empty-state-sub",
                        style={"marginBottom": "10px"},
                    ),
                    html.Label("Reference databases", className="form-label"),
                    html.Div(id="setup-db-recommendation", className="badge badge-blue", style={"marginBottom": "8px"}),
                    dcc.Checklist(
                        id="setup-db-sources",
                        options=[
                            {"label": " SILVA — rRNA (18S/16S)", "value": "silva"},
                            {"label": " PR2 — protists & eukaryotes", "value": "pr2"},
                            {"label": " BOLD — animals (COI barcodes)", "value": "bold"},
                            {"label": " NCBI — GenBank nucleotide", "value": "ncbi"},
                        ],
                        value=["silva", "pr2", "bold"],
                        className="checkbox-group source-checklist",
                        labelStyle={"display": "flex", "alignItems": "center", "gap": "4px", "marginBottom": "6px"},
                    ),
                    html.Hr(),
                    # ── Custom database upload (BLAST comparison) ──────────────
                    html.Label("Your own database (optional)", className="form-label"),
                    html.Div(
                        "Upload a FASTA file of reference sequences. Each of your OTUs/ASVs is compared "
                        "against it — matches show the species, and sequences with no match are flagged in red.",
                        className="empty-state-sub",
                        style={"marginBottom": "8px"},
                    ),
                    dcc.Upload(
                        id="setup-upload-custom-db",
                        children=html.Div(
                            [
                                html.Div("⬆ Drag & drop or click to choose a FASTA file", style={"fontWeight": "600"}),
                                html.Div(".fasta · .fa · .fasta.gz — up to 500 MB", className="empty-state-sub"),
                            ]
                        ),
                        className="upload-dropzone",
                        max_size=500 * 1024 * 1024,
                        multiple=False,
                    ),
                    dcc.Loading(html.Div(id="setup-custom-db-status", style={"marginTop": "8px"}), type="dot"),
                    html.Div(
                        id="db-params-container",
                        children=[
                            _advanced([
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
                                        value=_RUN_DEFAULTS.bold_params["mode"],
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
                            ]),
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
                    html.Span("5 · Build the tree", className="card-title"),
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
                    _advanced([
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
                        dcc.Input(id="setup-raxml-bootstraps", type="number", min=10, max=1000, step=10, value=_RUN_DEFAULTS.phylo_params["bs_trees"], className="dash-input"),
                        dcc.Checklist(
                            id="setup-adjust-direction",
                            options=[{"label": "Adjust direction", "value": "on"}],
                            value=["on"],
                            className="checkbox-group",
                        ),
                    ]),
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
            html.Div(id="setup-run-target", style={"marginBottom": "8px"}),
            html.Button("Run Pipeline", id="btn-run-pipeline", className="btn-run"),
        ],
        className="card card-body",
    )


def layout():
    return html.Div(
        [
            html.Div(
                [
                    html.Div("Set up your run", className="page-title"),
                    html.Div("Pick the steps to run, add your own database if you have one, then press Run.", className="page-sub"),
                ],
                className="page-header",
            ),
            html.Div(id="setup-active-project"),
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
                            html.Div(id="setup-preflight"),
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
