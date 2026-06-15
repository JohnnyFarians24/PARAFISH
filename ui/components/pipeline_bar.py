"""Horizontal pipeline progress bar."""

from __future__ import annotations

from dash import html


PIPELINE_STEPS = [
	("QC", "Quality control", "◫"),
	("Preprocess", "Trimmomatic + PEAR", "⟐"),
	("Analysis", "VSEARCH + DADA2", "⌬"),
	("Databases", "Taxonomy consensus", "◎"),
	("Phylo", "Trees and support", "☰"),
]


def pipeline_bar(step_states: list[str] | None = None):
	step_states = step_states or ["pending"] * len(PIPELINE_STEPS)
	steps = []
	for idx, ((name, desc, icon), state) in enumerate(zip(PIPELINE_STEPS, step_states), start=1):
		state_classes = []
		if state == "done":
			state_classes.append("is-done")
		if state == "active":
			state_classes.append("is-active")
		if state == "running":
			state_classes.append("is-running")
		if state == "error":
			state_classes.append("is-error")
		css = "pipe-step " + " ".join(state_classes) if state_classes else "pipe-step"
		steps.append(
			html.Div(
				[
					html.Span(str(idx), className="pipe-num"),
					html.Div(
						[
							html.Div(name, className="pipeline-step-label"),
							html.Div(desc, className="pipeline-step-subtitle"),
						],
						className="pipeline-step-copy",
					),
					html.Span(icon, className="sidebar-icon"),
				],
					className=css,
			)
		)
	return html.Div(steps, className="pipeline-bar")