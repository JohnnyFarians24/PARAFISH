"""Preprocess page.

Layer 1 (always visible): the KPI cards, the read-tracking table and the read-fate
plot. Read tracking is the table that tells you whether the run is trustworthy at all
— if a step silently ate 80% of the reads, that is the first thing anyone needs to
see, so it never goes behind a disclosure.

Layer 2 (collapsed): the parameter dump and the pipeline log. Detail you go looking
for, not detail you need to trust the numbers.
"""

from __future__ import annotations

from dash import html

from ui.components.collapsible import details
from ui.components.log_viewer import log_viewer


def layout():
	return html.Div(
		[
			html.Div([html.Div("Preprocess", className="page-title"), html.Div("Trimmomatic and PEAR parameters with per-sample progress.", className="page-sub")], className="page-header"),
			html.Div(id="preprocess-run-banner"),
			html.Div([html.Button("Run Preprocess", id="btn-run-preprocess", className="btn-run")], className="control-row"),
			html.Div(id="preprocess-summary", className="empty-state card"),
			html.Div(id="preprocess-pipeline-bar"),

			# ── Layer 1: did this work? ──────────────────────────────
			html.Div(
				[
					html.Div(id="preprocess-stat-input", className="stat-card"),
					html.Div(id="preprocess-stat-surviving", className="stat-card"),
					html.Div(id="preprocess-stat-assembly", className="stat-card"),
				],
				className="grid-3",
			),
			html.Div(
				[
					html.Div("Read tracking", className="card-title"),
					html.Div(
						"How many reads survive each step. Outside the initial quality filter no step should lose "
						"the majority of the reads — a red cell is a red flag. Paired-end counts are read pairs.",
						className="mini-note",
					),
					html.Div(id="preprocess-tracking-table"),
				],
				className="card card-body",
			),
			html.Div(
				[
					html.Div("Read fate per sample", className="card-title"),
					html.Div(id="preprocess-plots"),
				],
				className="card card-body",
			),

			# ── Layer 2: audit detail ────────────────────────────────
			details(
				"Parameters used",
				html.Div(
					[
						html.Div(id="preprocess-params"),
						html.Div("Adjust them in “Set up your run” before running the pipeline.", className="empty-state-sub", style={"marginTop": "8px"}),
					]
				),
				subtitle="Trimmomatic and PEAR",
			),
			details("Pipeline log", html.Div(log_viewer([]), id="preprocess-log-viewer")),
		],
		className="page-stack",
	)
