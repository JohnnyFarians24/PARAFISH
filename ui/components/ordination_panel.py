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
                    html.Label("Distância", className="mini-note"),
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
                    html.Label(f"Prevalência mínima (nº de amostras)", className="mini-note"),
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
                "As contagens de metabarcoding são composicionais: o total é imposto pela profundidade de "
                "sequenciação, não pela biologia. Por isso o método por omissão é a PCA de Aitchison (CLR) — "
                "um PCA sobre contagens brutas faria o PC1 captar a profundidade de sequenciação em vez de "
                "biologia.",
                className="mini-note",
            ),
            html.Div(
                f"Zeros substituídos por um pseudocount de {ZERO_REPLACEMENT} sobre as contagens. Nota: isto "
                f"torna o CLR apenas aproximadamente invariante à escala (o desvio é <1% mesmo para diferenças "
                f"de profundidade de 50x, muito abaixo do sinal biológico). A substituição multiplicativa "
                f"bayesiana (zCompositions) elimina esse resíduo, ao custo de uma dependência em R.",
                className="mini-note",
                style={"opacity": 0.8},
            ),
            controls(prefix, default_metric=default_metric),
            dcc.Checklist(
                id=f"{prefix}-exclude-host",
                options=[{"label": " Excluir hospedeiro e não classificados", "value": "exclude"}],
                value=["exclude"],
                className="checkbox-group",
            ),
            html.Div(
                "O sinal do hospedeiro domina a matriz e esmaga tudo o resto; sem esta exclusão o PC1 "
                "captura quase só a variação do hospedeiro.",
                className="mini-note",
            ),
            dcc.Loading(html.Div(id=f"{prefix}-plot"), type="dot"),
            html.Div(
                [
                    html.Button("Coordenadas (TSV)", id=f"{prefix}-download-coords", className="btn-parafish"),
                    html.Button("Matriz de distâncias (TSV)", id=f"{prefix}-download-dist", className="btn-parafish", style={"marginLeft": "8px"}),
                    dcc.Download(id=f"{prefix}-download-coords-data"),
                    dcc.Download(id=f"{prefix}-download-dist-data"),
                ],
                className="control-row",
            ),
            details("Biplot — que features causam a separação", dcc.Loading(html.Div(id=f"{prefix}-biplot"), type="dot")),
            details("Scree plot — variância por eixo", dcc.Loading(html.Div(id=f"{prefix}-scree"), type="dot")),
            details("PERMANOVA", dcc.Loading(html.Div(id=f"{prefix}-permanova"), type="dot")),
        ],
        style={"display": "flex", "flexDirection": "column", "gap": "8px"},
    )


def guard_message(n_samples: int, message: str):
    """The explicit empty state. Never a plot when the plot would mislead."""
    return html.Div(
        [
            html.Div("Ordenação indisponível", style={"fontWeight": "600", "color": "#d29922", "marginBottom": "6px"}),
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
        f"Atenção: {n_samples} amostras. A ordenação desenha-se, mas o poder interpretativo é muito "
        f"limitado abaixo de {LOW_CONFIDENCE_SAMPLES} amostras — trata a separação como indicativa, não como resultado.",
        className="mini-note",
        style={
            "padding": "8px 10px",
            "borderRadius": "8px",
            "background": "rgba(210,153,34,0.10)",
            "border": "1px solid rgba(210,153,34,0.35)",
            "marginBottom": "8px",
        },
    )
