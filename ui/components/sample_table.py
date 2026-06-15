"""Reusable sample table component."""

from __future__ import annotations

from dash import dcc, html


def _mode_badge(mode: str):
	mode_upper = str(mode).upper()
	css_class = "badge-teal" if mode_upper == "SE" else "badge-purple"
	return html.Span(mode_upper, className=f"badge {css_class}")


def _marker_badge(marker: str):
	marker_upper = str(marker).upper()
	css_map = {
		"COI": "badge-sand",
		"18S": "badge-amber",
		"ITS": "badge-blue",
		"28S": "badge-purple",
	}
	css_class = css_map.get(marker_upper, "badge-purple")
	return html.Span(marker_upper, className=f"badge {css_class}")


def render_sample_table(samples: list[dict] | None, selectable: bool = False):
	rows = samples or []
	headers = ["", "ID", "Label", "Mode", "Marker", "Status"] if selectable else ["ID", "Label", "Mode", "Marker", "Status"]
	head = html.Thead(html.Tr([html.Th(col) for col in headers]))
	body_rows = []
	for idx, sample in enumerate(rows):
		cells = []
		if selectable:
			cells.append(
				html.Td(
					dcc.Checklist(
						options=[{"label": "", "value": "selected"}],
						value=["selected"] if sample.get("selected") else [],
						id={"type": "sample-checkbox", "index": sample.get("internal_id") or sample.get("internal_sample_id") or str(idx)},
						style={"margin": "0"},
					)
				)
			)
		cells.extend(
			[
				html.Td(sample.get("internal_id", ""), className="table-id"),
				html.Td(sample.get("display_label", "")),
				html.Td(_mode_badge(sample.get("mode", "SE"))),
				html.Td(_marker_badge(sample.get("marker", "OTHER"))),
				html.Td(html.Span(sample.get("status", "pending"), className=f"badge badge-{sample.get('status_color', 'purple')}")),
			]
		)
		body_rows.append(html.Tr(cells))
	if not body_rows:
		body_rows.append(html.Tr(html.Td("No samples loaded yet.", colSpan=len(headers), className="empty-state-sub")))
	return html.Div(html.Table([head, html.Tbody(body_rows)], className="parafish-table"), className="card")