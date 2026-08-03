"""The ordination panel, parameterised so three pages can share one implementation.

Analysis owns the full controls (it is where the abundance matrix lives, and beta
diversity belongs next to alpha). Results shows a read-only compact version. Phylo
mounts the same panel with a UniFrac metric pre-selected. None of them re-implement it.
"""

from __future__ import annotations

from dash import dcc, html

from analysis.ordination import DISTANCES, LOW_CONFIDENCE_SAMPLES, MIN_PREVALENCE, MIN_SAMPLES, ZERO_REPLACEMENT
from ui.components.collapsible import details


def controls(prefix: str, *, default_metric: str = "aitchison"):
    """The metric / prevalence / host-exclusion controls."""
    return html.Div(
        [
            html.Div(
                [
                    html.Label("Distance", className="mini-note"),
                    dcc.Dropdown(
                        id=f"{prefix}-metric",
                        options=[{"label": label, "value": key} for key, label in DISTANCES.items()],
                        value=default_metric,
                        clearable=False,
                    ),
                ],
                style={"flex": "2", "minWidth": "260px"},
            ),
            html.Div(
                [
                    html.Label("Colorir por", className="mini-note"),
                    dcc.Dropdown(id=f"{prefix}-colour", options=[], value="group", clearable=False),
                ],
                style={"flex": "1", "minWidth": "160px"},
            ),
            html.Div(
                [
                    html.Label(f"Minimum prevalence (no. of samples)", className="mini-note"),
                    dcc.Slider(
                        id=f"{prefix}-prevalence",
                        min=1,
                        max=10,
                        step=1,
                        value=MIN_PREVALENCE,
                        marks={1: "1", 2: "2", 5: "5", 10: "10"},
                        tooltip={"placement": "bottom", "always_visible": False},
                    ),
                ],
                style={"flex": "1", "minWidth": "200px"},
            ),
        ],
        style={"display": "flex", "gap": "16px", "flexWrap": "wrap", "marginBottom": "10px"},
    )


def panel(prefix: str, *, compact: bool = False, default_metric: str = "aitchison"):
    """The full panel. ``compact=True`` gives Results the plot only, no controls."""
    if compact:
        return html.Div(
            [
                dcc.Loading(html.Div(id=f"{prefix}-plot"), type="dot"),
                html.Div(
                    html.A("Controlos completos, biplot e PERMANOVA em Analysis →", href="/analysis", className="mini-note"),
                    style={"marginTop": "6px"},
                ),
            ]
        )

    return html.Div(
        [
            html.Div(
                "Metabarcoding counts are compositional: the total is imposed by sequencing depth, "
                "not by biology. That is why the default method is Aitchison PCA (CLR) — "
                "a PCA on raw counts would make PC1 capture sequencing depth instead of "
                "biology.",
                className="mini-note",
            ),
            html.Div(
                f"Zeros replaced by a pseudocount of {ZERO_REPLACEMENT} on the counts. Note: this "
                f"makes the CLR only approximately scale-invariant (the drift is <1% even for depth "
                f"differences of 50x, far below the biological signal). Bayesian multiplicative "
                f"replacement (zCompositions) removes that residue, at the cost of an R dependency.",
                className="mini-note",
                style={"opacity": 0.8},
            ),
            controls(prefix, default_metric=default_metric),
            dcc.Checklist(
                id=f"{prefix}-exclude-host",
                options=[{"label": " Exclude host and unclassified", "value": "exclude"}],
                value=["exclude"],
                className="checkbox-group",
            ),
            html.Div(
                "The host signal dominates the matrix and swamps everything else; without this exclusion PC1 "
                "captures almost nothing but host variation.",
                className="mini-note",
            ),
            dcc.Loading(html.Div(id=f"{prefix}-plot"), type="dot"),
            html.Div(
                [
                    html.Button("Coordenadas (TSV)", id=f"{prefix}-download-coords", className="btn-parafish"),
                    html.Button("Distance matrix (TSV)", id=f"{prefix}-download-dist", className="btn-parafish", style={"marginLeft": "8px"}),
                    dcc.Download(id=f"{prefix}-download-coords-data"),
                    dcc.Download(id=f"{prefix}-download-dist-data"),
                ],
                className="control-row",
            ),
            details("Biplot — which features drive the separation", dcc.Loading(html.Div(id=f"{prefix}-biplot"), type="dot")),
            details("Scree plot — variance per axis", dcc.Loading(html.Div(id=f"{prefix}-scree"), type="dot")),
            details("PERMANOVA", dcc.Loading(html.Div(id=f"{prefix}-permanova"), type="dot")),
        ],
        style={"display": "flex", "flexDirection": "column", "gap": "8px"},
    )


def guard_message(n_samples: int, message: str):
    """The explicit empty state. Never a plot when the plot would mislead."""
    return html.Div(
        [
            html.Div("Ordination unavailable", style={"fontWeight": "600", "color": "#d29922", "marginBottom": "6px"}),
            html.Div(message, className="empty-state-sub"),
        ],
        style={
            "padding": "14px 16px",
            "borderRadius": "8px",
            "background": "rgba(210,153,34,0.08)",
            "border": "1px solid rgba(210,153,34,0.30)",
        },
    )


def low_confidence_warning(n_samples: int):
    """Rendered above the plot when it is drawable but barely interpretable."""
    return html.Div(
        f"Warning: {n_samples} samples. The ordination will render, but interpretive power is very "
        f"limited below {LOW_CONFIDENCE_SAMPLES} samples — read the separation as indicative, not as a result.",
        className="mini-note",
        style={
            "padding": "8px 10px",
            "borderRadius": "8px",
            "background": "rgba(210,153,34,0.10)",
            "border": "1px solid rgba(210,153,34,0.35)",
            "marginBottom": "8px",
        },
    )
