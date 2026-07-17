"""Analysis page."""

from __future__ import annotations

from dash import dcc, html

from ui.components.collapsible import details
from ui.components.ordination_panel import panel as ordination_panel


def layout():
	return html.Div(
		[
			html.Div([html.Div("Analysis", className="page-title"), html.Div("VSEARCH and DADA2 side by side.", className="page-sub")], className="page-header"),
			html.Div([html.Button("Run Analysis", id="btn-run-analysis", className="btn-run")], className="control-row"),
			html.Div(
				[
					html.Div(id="analysis-summary", className="empty-state card"),
					html.Div(id="analysis-pipeline-bar"),
					html.Div(
						[
							html.Div(id="analysis-stat-otus", className="stat-card"),
							html.Div(id="analysis-stat-asvs", className="stat-card"),
							html.Div(id="analysis-stat-consensus", className="stat-card"),
						],
						className="grid-3",
					),
				],
				className="card card-body",
			),
			# ── Layer 1: the numbers that establish trust ────────────
			html.Div(
				[
					html.Div("Per-sample yield", className="card-title"),
					html.Div("Reads surviving filtering and the OTUs/ASVs recovered from each sample.", className="mini-note"),
					html.Div(id="analysis-sample-table"),
				],
				className="card card-body",
			),
			html.Div(
				[
					html.Div("Alpha diversity", className="card-title"),
					html.Div(id="analysis-rarefaction-note", className="mini-note"),
					html.Div(id="analysis-alpha-table"),
				],
				className="card card-body",
			),
			# Beta sits directly under alpha: within-sample and between-sample diversity
			# are the two halves of the same analysis, and this is the page that owns the
			# abundance matrix both are computed from.
			html.Div(
				[
					html.Div("Beta diversity — ordination", className="card-title"),
					ordination_panel(prefix="ord-analysis"),
				],
				className="card card-body",
			),

			# ── Layer 2: audit detail ────────────────────────────────
			details(
				"Chimeras",
				html.Div(
					[
						html.Div(
							"Reported two ways on purpose: a run can flag many chimeric variants that together "
							"carry almost no reads. Neither percentage implies the other.",
							className="mini-note",
						),
						html.Div(id="analysis-chimera"),
					]
				),
			),
			details("DADA2 — read tracking", html.Div(id="analysis-dada2-tracking")),
			html.Div(
				[
					# Tabs, not an accordion: absolute and relative are two mutually
					# exclusive views of the SAME matrix, which is exactly what tabs are for.
					html.Div("Abundance table (OTU / ASV x sample)", className="card-title"),
					dcc.RadioItems(
						id="analysis-abundance-mode",
						options=[
							{"label": " Absolute counts", "value": "absolute"},
							{"label": " Relative abundance (%)", "value": "relative"},
						],
						value="absolute",
						inline=True,
						className="checkbox-group",
					),
					html.Div(
						"Compare composition between samples with relative abundance; filter rare features on "
						"absolute counts.",
						className="mini-note",
					),
					html.Div(id="analysis-otu-table"),
				],
				className="card card-body",
			),
			# Layer 1 hero: rank-abundance is the one picture that shows whether the
			# community is dominated by a handful of features or evenly spread.
			html.Div(
				[
					html.Div("Diversity and composition", className="card-title"),
					dcc.Loading(html.Div(id="analysis-plots"), type="dot"),
				],
				className="card card-body",
			),
			details("Shared and exclusive features", dcc.Loading(html.Div(id="analysis-upset"), type="dot")),
			html.Div(
				[
					html.Div("Which samples each OTU/ASV appears in", className="card-title"),
					html.Div(
						"One row per feature: in how many and in WHICH samples it was found. Those shared "
						"by several samples come first. Each set (OTU/ASV) has its own table.",
						className="mini-note",
					),
					dcc.RadioItems(
						id="analysis-occupancy-set",
						options=[{"label": " OTUs", "value": "otu"}, {"label": " ASVs", "value": "asv"}],
						value="otu",
						inline=True,
						className="chip-list",
						inputStyle={"marginRight": "6px"},
						labelStyle={"marginRight": "18px"},
					),
					dcc.RadioItems(
						id="analysis-occupancy-filter",
						options=[
							{"label": " All", "value": "all"},
							{"label": " Shared only (>1 sample)", "value": "shared"},
							{"label": " Exclusive only (1 sample)", "value": "unique"},
						],
						value="shared",
						inline=True,
						className="chip-list",
						inputStyle={"marginRight": "6px"},
						labelStyle={"marginRight": "18px"},
					),
					dcc.Loading(html.Div(id="analysis-occupancy"), type="dot"),
				],
				className="card card-body",
			),
			html.Div(
				[
					html.Div("Downloads", className="card-title"),
					html.Div(id="analysis-downloads"),
				],
				className="card card-body",
			),
		],
		className="page-stack",
	)
