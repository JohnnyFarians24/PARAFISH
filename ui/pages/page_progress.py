"""Pipeline progress page."""

from __future__ import annotations

from dash import dcc, html


def layout():
    return html.Div(
        [
            html.Div("Pipeline Progress", className="page-title"),
            html.Div("Live execution and logs", className="page-sub"),
            dcc.Interval(id="progress-interval", interval=2000, n_intervals=0, disabled=True),
            html.Div(id="progress-overall-status", className="badge badge-blue"),
            html.Button("Start Pipeline", id="btn-start-progress", className="btn-run"),
            html.Div(id="progress-stage-rows"),
            html.Div(id="progress-logs", className="monospace"),
        ],
        className="main-content",
    )
