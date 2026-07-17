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
	base_headers = ["ID", "Label", "Mode", "Marker", "Status"]
	# Negative controls and experimental groups are only editable in the selectable
	# view (where the run is configured); the read-only view just lists the samples.
	headers = ["", *base_headers, "Control", "Group"] if selectable else base_headers
	head = html.Thead(html.Tr([html.Th(col) for col in headers]))
	body_rows = []
	for idx, sample in enumerate(rows):
		sample_id = sample.get("internal_id") or sample.get("internal_sample_id") or str(idx)
		cells = []
		if selectable:
			cells.append(
				html.Td(
					dcc.Checklist(
						options=[{"label": "", "value": "selected"}],
						value=["selected"] if sample.get("selected") else [],
						id={"type": "sample-checkbox", "index": sample_id},
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
		if selectable:
			cells.append(
				html.Td(
					dcc.Checklist(
						options=[{"label": "", "value": "control"}],
						value=["control"] if sample.get("is_control") else [],
						id={"type": "sample-control", "index": sample_id},
						style={"margin": "0"},
					),
					title="Negative / blank control: excluded from diversity, used to call contaminants",
				)
			)
			cells.append(
				html.Td(
					dcc.Input(
						id={"type": "sample-group", "index": sample_id},
						value=sample.get("group") or "",
						placeholder="e.g. site A",
						debounce=True,
						style={"width": "110px"},
					),
					title="Experimental group (site, host, treatment) used for grouped comparisons",
				)
			)
		body_rows.append(html.Tr(cells))
	if not body_rows:
		body_rows.append(html.Tr(html.Td("No samples loaded yet.", colSpan=len(headers), className="empty-state-sub")))
	return html.Div(html.Table([head, html.Tbody(body_rows)], className="parafish-table"), className="card")