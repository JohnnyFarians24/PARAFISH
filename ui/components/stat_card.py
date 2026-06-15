"""Stat card component."""

from __future__ import annotations

from dash import html


def stat_card(value: str, label: str, color: str = "var(--text-0)"):
	return html.Div(
		[
			html.Div(value, className="stat-value", style={"color": color}),
			html.Div(label, className="stat-label"),
		],
		className="stat-card",
	)