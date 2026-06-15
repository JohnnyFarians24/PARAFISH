"""Analysis page."""

from __future__ import annotations

from dash import html


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
		],
		className="page-stack",
	)