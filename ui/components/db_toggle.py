"""Database activation toggle component."""

from __future__ import annotations

from dash import dcc, html


def db_toggle(label: str, value: bool = False, toggle_id: str | None = None):
	effective_id = toggle_id or f"db-toggle-{label.lower().replace(' ', '-')}"
	return html.Div(
		[
			html.Div(label, className="db-toggle-name"),
			html.Div("active" if value else "inactive", className="db-toggle-version"),
			dcc.Checklist(
				id=effective_id,
				options=[{"label": label, "value": "enabled"}],
				value=["enabled"] if value else [],
				className="db-toggle",
			),
		],
		className="db-toggle-row",
	)