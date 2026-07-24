"""Phylogenetics page.

Layer 1: stage status, KPI cards, and the annotated tree — the tree IS the result here,
so it is never behind a click.

Layer 2: alignment QC, branch support, UniFrac ordination, phylogenetic alpha diversity.
"""

from __future__ import annotations

from dash import dcc, html

from ui.components.collapsible import details
from ui.components.ordination_panel import panel as ordination_panel


def layout():
	return html.Div(
		[
			html.Div(
				[
					html.Div("Phylogenetics", className="page-title"),
					html.Div("FastTree, RAxML-NG and hybrid routing for aligned OTUs/ASVs.", className="page-sub"),
				],
				className="page-header",
			),
			html.Div(id="phylo-run-banner"),
			html.Div([html.Button("Run Phylo", id="btn-run-phylo", className="btn-run")], className="control-row"),
			html.Div(id="phylo-summary", className="empty-state card"),
			html.Div(id="phylo-pipeline-bar"),

			# ── Layer 1 ──────────────────────────────────────────────
			html.Div(
				[
					html.Div(id="phylo-stat-tips", className="stat-card"),
					html.Div(id="phylo-stat-alignment", className="stat-card"),
					html.Div(id="phylo-stat-support", className="stat-card"),
				],
				className="grid-3",
			),
			html.Div(
				[
					html.Div("Annotated tree", className="card-title"),
					html.Div(id="phylo-tree-viewer"),
				],
				className="card card-body",
			),

			# ── Layer 2 ──────────────────────────────────────────────
			details(
				"Phylogenetic ordination (UniFrac)",
				html.Div(
					[
						html.Div(
							"Complements Bray-Curtis: UniFrac measures how much evolutionary history two samples "
							"share, not just how many features. Weighted uses abundance; unweighted uses "
							"presence/absence.",
							className="mini-note",
						),
						# Same reusable panel as Analysis, with a phylogenetic metric preselected.
						ordination_panel(prefix="ord-phylo", default_metric="weighted_unifrac"),
					]
				),
			),
			details(
				"Alignment QC (MAFFT)",
				html.Div(
					[
						html.Div(
							"Conservation and gaps together: a column that is 90% gaps but whose single "
							"remaining base is invariant scores perfectly if you look at conservation alone.",
							className="mini-note",
						),
						dcc.Loading(html.Div(id="phylo-alignment"), type="dot"),
					]
				),
			),
			details("Branch support", dcc.Loading(html.Div(id="phylo-support"), type="dot")),
			details("Phylogenetic alpha diversity", html.Div(id="phylo-visuals")),
			# Retained for the existing callback; superseded by the ordination panel above.
			html.Div(id="phylo-unifrac", style={"display": "none"}),
		],
		className="page-stack",
	)
