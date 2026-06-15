"""Progress bar component."""

from __future__ import annotations

from dash import html


def progress_bar(value: float, label: str = ""):
	value = max(0.0, min(100.0, float(value)))
	return html.Div(
		[
			html.Div(label, className="muted-text"),
			html.Div(
				html.Div(style={"width": f"{value}%", "height": "100%", "background": "var(--blue)", "borderRadius": "999px"}),
				style={"height": "8px", "background": "var(--bg-4)", "borderRadius": "999px", "overflow": "hidden"},
			),
		],
		style={"display": "flex", "flexDirection": "column", "gap": "6px"},
	)