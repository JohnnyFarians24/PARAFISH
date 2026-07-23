"""Pipeline progress page."""

from __future__ import annotations

from dash import dcc, html


def layout():
    return html.Div(
        [
            html.Div("Pipeline Progress", className="page-title"),
            html.Div("Configure on Setup, then press Start Pipeline to run", className="page-sub"),
            dcc.Interval(id="progress-interval", interval=2000, n_intervals=0, disabled=False),
            # Idle text, not a blank badge: reaching this page no longer means a run has
            # started (Setup only navigates here now), so the page must say it is waiting
            # for the click instead of looking like a frozen, blocked screen.
            html.Div(
                "Ready — press Start Pipeline to begin the run.",
                id="progress-overall-status",
                className="badge badge-blue",
            ),
            html.Button("Start Pipeline", id="btn-start-progress", className="btn-run"),
            html.Div(id="progress-stage-rows"),
            html.Div(id="progress-logs", className="monospace"),
        ],
        className="main-content",
    )
