"""Status badge component."""

from __future__ import annotations

from dash import html


STATUS_CLASSES = {
	"done": "badge-green",
	"success": "badge-green",
	"active": "badge-blue",
	"running": "badge-amber",
	"pending": "badge-purple",
	"error": "badge-red",
	"failure": "badge-red",
	"blocked": "badge-red",
}


def status_badge(label: str, status: str = "pending"):
	css_class = STATUS_CLASSES.get(status, "badge-purple")
	return html.Span(label, className=f"badge {css_class}")