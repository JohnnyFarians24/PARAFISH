"""QC page.

Layer 1 (always visible): the KPI cards, the run-overview table and the read-depth
plot. Those three answer "is this data usable?" without a click.

Layer 2 (collapsed): the per-base quality curves and the length distribution — you go
to them when the headline says something is off.
"""

from __future__ import annotations

from dash import html

from ui.components.collapsible import details


def layout():
	return html.Div(
		[
			html.Div([html.Div("Quality control", className="page-title"), html.Div("FastQC summaries before trimming.", className="page-sub")], className="page-header"),
			html.Div(id="qc-run-banner"),
			html.Div([html.Button("Run QC", id="btn-run-qc", className="btn-run")], className="control-row"),
			html.Div(id="qc-summary", className="empty-state card"),

			# ── Layer 1 ──────────────────────────────────────────────
			html.Div(
				[
					html.Div(id="qc-stat-total", className="stat-card"),
					html.Div(id="qc-stat-pass", className="stat-card"),
					html.Div(id="qc-stat-reads", className="stat-card"),
					html.Div(id="qc-stat-gc", className="stat-card"),
				],
				className="grid-4",
			),
			html.Div(
				[
					html.Div(
						[
							html.Div("Run overview", className="card-title"),
							html.Div(id="qc-sample-table"),
						],
						className="card-body",
					),
				],
				className="card",
			),
			html.Div(
				[
					html.Div("Read depth", className="card-title"),
					html.Div(id="qc-plot-hero"),
				],
				className="card card-body",
			),

			# ── Layer 2 ──────────────────────────────────────────────
			details(
				"Quality curves and length distribution",
				html.Div(id="qc-plots"),
				subtitle="per-position quality, read lengths",
			),
		],
		className="page-stack",
	)