"""Results — the executive summary of a run.

Results is not a dump of what the other pages already show. It answers three questions,
in this order, and a reader should never need to leave it to answer them:

1. *Is this run valid?*   — run header, read funnel, per-stage KPI cards.
2. *What was found?*      — community composition, ordination.
3. *Where are the files?* — exports, with anything missing explained rather than hidden.

Layer 1 is everything above the fold. Layer 2 is the full tables. Layer 3 is the log and
the reproducibility block.
"""

from __future__ import annotations

from dash import dcc, html

from ui.components.collapsible import details
from ui.components.ordination_panel import panel as ordination_panel


def _stage_viewer():
    """The per-stage / per-sample browser, kept but demoted below the summary."""
    return details(
        "Explore outputs by stage and sample",
        html.Div(
            [
                html.Div(
                    [
                        dcc.Dropdown(
                            id="results-stage-selector",
                            placeholder="Choose a stage...",
                            style={"width": "260px", "display": "inline-block", "marginRight": "8px"},
                        ),
                        dcc.Dropdown(
                            id="results-sample-selector",
                            placeholder="Choose a sample...",
                            style={"width": "260px", "display": "inline-block"},
                            searchable=True,
                        ),
                    ],
                    style={"marginBottom": "12px"},
                ),
                html.Div(id="results-viewer"),
            ]
        ),
    )


def layout():
    return html.Div(
        [
            html.Div(
                [
                    html.Div("Results", className="page-title"),
                    html.Div("Run summary: validity, what was found, and where the files live.", className="page-sub"),
                ],
                className="page-header",
            ),

            # ── Bloco A — Run header (Layer 1) ───────────────────────
            html.Div(id="results-run-header", className="card card-body"),

            # ── Bloco B — Read funnel (Layer 1, hero) ────────────────
            html.Div(
                [
                    html.Div("End-to-end read funnel", className="card-title"),
                    html.Div(
                        "Every read, from raw FASTQ to classified feature. A large drop at a single "
                        "step is the first place to investigate — steps that lose more than 20% are shown in red "
                        "(except the initial quality filter, where losing reads is the point).",
                        className="mini-note",
                    ),
                    dcc.Loading(html.Div(id="results-funnel"), type="dot"),
                ],
                className="card card-body",
            ),

            # ── Bloco C — Per-stage KPI cards (Layer 1) ──────────────
            html.Div(
                [
                    html.Div("Per-stage summary", className="card-title"),
                    html.Div("Each card links to its page.", className="mini-note"),
                    dcc.Loading(html.Div(id="results-stage-cards"), type="dot"),
                ],
                className="card card-body",
            ),

            # ── Bloco D — Community composition (Layer 1) ────────────
            html.Div(
                [
                    html.Div("Community composition", className="card-title"),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Label("Taxonomic level", className="mini-note"),
                                    dcc.Dropdown(
                                        id="results-composition-rank",
                                        options=[
                                            {"label": "Kingdom", "value": "Consensus_Kingdom"},
                                            {"label": "Phylum", "value": "Consensus_Phylum"},
                                            {"label": "Class", "value": "Consensus_Class"},
                                            {"label": "Order", "value": "Consensus_Order"},
                                            {"label": "Family", "value": "Consensus_Family"},
                                            {"label": "Genus", "value": "Consensus_Genus"},
                                        ],
                                        value="Consensus_Phylum",
                                        clearable=False,
                                    ),
                                ],
                                style={"flex": "1", "minWidth": "200px"},
                            ),
                            html.Div(
                                dcc.Checklist(
                                    id="results-composition-exclude",
                                    options=[{"label": " Exclude host and unclassified", "value": "exclude"}],
                                    value=[],
                                    className="checkbox-group",
                                ),
                                style={"flex": "1", "minWidth": "260px", "alignSelf": "end"},
                            ),
                        ],
                        style={"display": "flex", "gap": "16px", "flexWrap": "wrap", "marginBottom": "10px"},
                    ),
                    dcc.Loading(html.Div(id="results-composition"), type="dot"),
                ],
                className="card card-body",
            ),

            # ── Bloco E — Ordination, compact (Layer 1) ──────────────
            html.Div(
                [
                    html.Div("Ordination (beta diversity)", className="card-title"),
                    ordination_panel(prefix="ord-results", compact=True),
                ],
                className="card card-body",
            ),

            # ── Bloco F — Tables (Layer 2) ───────────────────────────
            details("Alpha diversity per sample", dcc.Loading(html.Div(id="results-alpha-table"), type="dot")),
            details("Taxonomy table (per-feature consensus)", dcc.Loading(html.Div(id="results-taxonomy-table"), type="dot")),
            details("Detailed read tracking", dcc.Loading(html.Div(id="results-tracking-table"), type="dot")),
            details("Between-sample comparison", html.Div(id="results-comparison")),
            details(
                "Target taxa (parasites)",
                html.Div(
                    [
                        dcc.Input(
                            id="results-target-input",
                            type="text",
                            placeholder="e.g. Diplostomum, Diplostomum spathaceum",
                            debounce=True,
                            style={"width": "100%", "maxWidth": "480px"},
                        ),
                        html.Div(
                            "Comma-separated taxa: where each one appears, its prevalence with 95% CI, "
                            "and how well supported the assignment is.",
                            className="mini-note",
                        ),
                        dcc.Loading(html.Div(id="results-targets"), type="dot"),
                    ]
                ),
            ),
            _stage_viewer(),

            # ── Bloco G — Exports (Layer 2) ──────────────────────────
            details(
                "Downloads",
                html.Div(
                    [
                        html.Button("Download all (zip)", id="btn-download-all", className="btn-parafish"),
                        dcc.Download(id="results-download-all"),
                        dcc.Loading(html.Div(id="results-exports"), type="dot"),
                    ]
                ),
                open_by_default=True,
            ),

            # ── Bloco H — Reproducibility (Layer 3) ──────────────────
            details("Reproducibility — tools, parameters and versions", dcc.Loading(html.Div(id="results-repro"), type="dot")),
            details("Run log", html.Div(id="results-logs", className="monospace")),

            html.Div(html.Button("New Run", id="btn-new-run", className="btn-parafish"), style={"marginTop": "16px"}),
            # Kept for callbacks that still target them; the summary above supersedes them.
            html.Div(id="results-sankey", style={"display": "none"}),
            html.Div(id="results-dashboard", style={"display": "none"}),
            html.Div(id="results-stats", style={"display": "none"}),
        ],
        className="page-stack",
    )
