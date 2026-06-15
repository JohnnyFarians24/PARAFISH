"""Results page."""

from __future__ import annotations

from dash import dcc, html


def layout():
	return html.Div(
		[
			html.Div([html.Div("Results", className="page-title"), html.Div("Final run summary and export hub.", className="page-sub")], className="page-header"),
			html.Button("New Run", id="btn-new-run", className="btn-parafish", style={"marginTop": "16px"}),
			html.Div(
				[
					dcc.Dropdown(id="results-stage-selector", placeholder="Select stage...", style={"width": "200px", "display": "inline-block", "marginRight": "8px"}),
						dcc.Dropdown(id="results-sample-selector", placeholder="Select sample...", style={"width": "200px", "display": "inline-block"}, searchable=True),
				],
				style={"marginBottom": "12px"},
			),
			html.Div(id="results-viewer", className="card card-body"),
			html.Div(id="results-stats", className="grid-3 card-body"),
			html.Div(id="results-exports", className="card-body"),
			html.Div(id="results-logs", className="card-body monospace"),
		],
		className="page-stack",
	)