"""Phylogenetics page."""

from __future__ import annotations

from dash import html

def layout():
	return html.Div(
		[
			html.Div(
				[
					html.Div("Phylogenetics", className="page-title"),
					html.Div(
						"FastTree, RAxML-NG and hybrid routing for aligned ASVs.",
						className="page-sub",
					),
				],
				className="page-header",
			),
			html.Div(
				[
					html.Button("Run Phylo", id="btn-run-phylo", className="btn-run"),
				],
				className="control-row",
			),
			html.Div(
				[
					html.Div(id="phylo-summary", className="empty-state card"),
					html.Div(id="phylo-pipeline-bar"),
					html.Div(id="phylo-tree-viewer", className="empty-state-sub card-body"),
				],
				className="card",
			),
		],
		className="page-stack",
	)