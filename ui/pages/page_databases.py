"""Databases page."""

from __future__ import annotations

from dash import dcc, html

from ui.components.db_toggle import db_toggle


def layout():
	return html.Div(
		[
			html.Div([html.Div("Databases", className="page-title"), html.Div("Enable and inspect SILVA, PR2, NCBI and BOLD sources.", className="page-sub")], className="page-header"),
			html.Div([html.Button("Run Databases", id="btn-run-databases", className="btn-run")], className="control-row"),
			html.Div(
				[
					html.Div([
						html.Div("Active databases", className="card-title"),
						html.Div(id="databases-summary", className="empty-state card"),
						db_toggle("SILVA 138.2 · SSU", True, "db-silva"),
						db_toggle("PR2 5.1.1", True, "db-pr2"),
						db_toggle("BOLD", True, "db-bold"),
						db_toggle("NCBI nt", False, "db-ncbi"),
						html.Hr(),
						html.Div("Custom Database", className="card-title"),
						dcc.Upload(
							id="upload-custom-db",
							children=html.Div(["Drag and drop or ", html.A("select a FASTA file"), " (.fasta, .fa, .fasta.gz, max 500 MB)"]),
							style={"borderStyle": "dashed", "borderWidth": "1px", "borderRadius": "4px", "padding": "12px", "textAlign": "center"},
							max_size=500 * 1024 * 1024,
							multiple=False,
						),
						html.Div(id="custom-db-upload-status"),
					], className="card card-body"),
					html.Div([html.Div("Consensus table and per-source tabs will appear here.", className="empty-state-sub")], className="card card-body"),
				],
				className="grid-2",
			),
		],
		className="page-stack",
	)