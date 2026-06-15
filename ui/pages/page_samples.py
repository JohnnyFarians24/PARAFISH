"""Samples page."""

from __future__ import annotations

from dash import dcc, html

from ui.components.sample_table import render_sample_table


def layout():
	return html.Div(
		[
			html.Div([html.Div("Samples", className="page-title"), html.Div("Upload FASTQ files and inspect the detected sample structure.", className="page-sub")], className="page-header"),
			html.Div(
				[
				dcc.Upload(
					id="upload-fastq",
					multiple=True,
					max_size=2 * 1024 * 1024 * 1024,
					accept=".fastq,.fastq.gz,.fq,.fq.gz,.zip,.tar.gz,.rar",
					children=html.Div(
						[
							html.Div("⇪", className="upload-zone-icon"),
							html.Div("Drag and drop FASTQ files, a ZIP archive, or browse", className="upload-zone-text"),
							html.Div("Accepted: .fastq, .fastq.gz, .fq, .fq.gz, .zip, .rar", className="upload-zone-sub"),
						],
						className="upload-zone",
					),
					className="upload-zone",
				),
				html.Div([html.Div("Ready for upload", className="badge badge-green"), html.Div("The sample detector will split PE runs and preserve SE routing.", className="empty-state-sub")], className="control-row"),
					html.Button("Clear samples", id="btn-clear-samples", className="btn-parafish", style={"marginTop": "8px"}),
				],
			className="card card-body",
			),
			html.Div(id="samples-upload-status", className="chip-list"),
			html.Div(id="samples-table", children=render_sample_table([], selectable=True)),
		],
		className="page-stack",
	)