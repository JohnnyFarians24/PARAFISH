"""Sidebar navigation component."""

from __future__ import annotations

from dash import dcc, html


STEP_STATES = {
	"done": "var(--green)",
	"active": "var(--blue)",
	"running": "var(--amber)",
	"pending": "var(--bg-4)",
	"error": "var(--red)",
}


NAV_ITEMS = [
	("/project", "Project", "⌂", "Workspace and session"),
	("/samples", "Samples", "⇪", "Upload and detection"),
	("/setup", "Setup", "⚙", "Pipeline configuration"),
	("/qc", "QC", "◫", "FastQC summaries"),
	("/preprocess", "Preprocess", "⟐", "Trimming and merging"),
	("/analysis", "Analysis", "⌬", "ASV inference and clustering"),
	("/databases", "Databases", "◎", "Taxonomy references"),
	("/phylo", "Phylo", "☰", "Tree building"),
	("/progress", "Progress", "▶", "Pipeline execution"),
	("/results", "Results", "⟡", "Exports and reports"),
]


def render_sidebar(active_path: str = "/project", step_states: list[str] | None = None):
	step_states = step_states or ["active"] + ["pending"] * (len(NAV_ITEMS) - 1)
	always_accessible = {"/project", "/samples", "/setup"}
	items = []
	for (href, label, icon, subtitle), state in zip(NAV_ITEMS, step_states):
		is_disabled = state == "pending" and href != active_path and href not in always_accessible
		items.append(
			dcc.Link(
				[
					html.Span(className="step-dot"),
					html.Span(icon, className="sidebar-icon"),
					html.Div(
						[
							html.Div(label, className="sidebar-label"),
							html.Div(subtitle, className="sidebar-subtitle"),
						],
						className="sidebar-copy",
					),
				],
				href=href,
				className=f"sidebar-item {'active' if href == active_path else ''} {'is-done' if state == 'done' else ''} {'is-running' if state == 'running' else ''} {'is-disabled' if is_disabled else ''} {'is-error' if state == 'error' else ''}",
			)
		)
	return html.Div(
		[
			html.Div("Workflow", className="sidebar-section"),
			*items,
			html.Div(
				[
					html.Div("Parafish shell", className="project-name-sidebar"),
					html.Div("A phased sequencing workspace", className="project-meta-sidebar"),
				],
				className="sidebar-footer",
			),
		],
	)