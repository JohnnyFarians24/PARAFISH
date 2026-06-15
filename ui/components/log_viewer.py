"""Log viewer component."""

from __future__ import annotations

from dash import html


LEVEL_CLASSES = {
	"OK": "badge-green",
	"RUN": "badge-amber",
	"INF": "badge-blue",
	"ERR": "badge-red",
}


def log_viewer(log_lines: list[dict] | None):
	log_lines = log_lines or []
	rows = []
	for line in log_lines[-50:]:
		level = str(line.get("level", "INF")).upper()
		rows.append(
			html.Div(
				[
					html.Span(line.get("ts", "--:--:--"), className="log-line-ts"),
					html.Span(level, className=f"badge {LEVEL_CLASSES.get(level, 'badge-purple')}"),
					html.Span(line.get("msg", ""), className=f"log-line-{level.lower()[:3] if level else 'inf'}"),
				],
				className="log-line",
			)
		)
	if not rows:
		rows = [html.Div("No log lines yet.", className="empty-state-sub")]
	return html.Div(rows, className="log-viewer")