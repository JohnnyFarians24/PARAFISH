"""Preprocess page."""

from __future__ import annotations

from dash import html

from ui.components.log_viewer import log_viewer


def layout():
	return html.Div(
		[
			html.Div([html.Div("Preprocess", className="page-title"), html.Div("Trimmomatic and PEAR parameters with per-sample progress.", className="page-sub")], className="page-header"),
			html.Div([html.Button("Run Preprocess", id="btn-run-preprocess", className="btn-run")], className="control-row"),
			html.Div(
				[
					html.Div([html.Div("Parameters", className="card-title"), html.Div(id="preprocess-summary", className="empty-state card"), html.Div("Configure Trimmomatic and PEAR here.", className="empty-state-sub")], className="card card-body"),
					html.Div([html.Div(id="preprocess-pipeline-bar"), html.Div(log_viewer([]), id="preprocess-log-viewer")], className="card card-body"),
				],
				className="grid-2",
			),
		],
		className="page-stack",
	)