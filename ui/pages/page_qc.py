"""QC page."""

from __future__ import annotations

from dash import html


def layout():
	return html.Div(
		[
			html.Div([html.Div("Quality control", className="page-title"), html.Div("FastQC summaries before trimming.", className="page-sub")], className="page-header"),
			html.Div([html.Button("Run QC", id="btn-run-qc", className="btn-run")], className="control-row"),
			html.Div(id="qc-summary", className="empty-state card"),
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
		],
		className="page-stack",
	)