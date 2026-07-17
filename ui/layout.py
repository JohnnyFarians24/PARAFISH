"""Dash layout definition for the modular scientific workflow UI."""

from __future__ import annotations

from importlib import import_module

from ui.components.connection_banner import connection_banner
from ui.components.sidebar import render_sidebar
from ui.components.topbar import render_topbar
from ui.pages.page_analysis import layout as analysis_layout
from ui.pages.page_databases import layout as databases_layout
from ui.pages.page_phylo import layout as phylo_layout
from ui.pages.page_preprocess import layout as preprocess_layout
from ui.pages.page_project import layout as project_layout
from ui.pages.page_qc import layout as qc_layout
from ui.pages.page_results import layout as results_layout
from ui.pages.page_samples import layout as samples_layout
from ui.pages.page_setup import layout as setup_layout
from ui.pages.page_progress import layout as progress_layout


def create_layout():
    """Return the full application shell with stores and page container."""
    try:
        dash_module = import_module("dash")
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Dash is required to build the UI layout.") from exc

    dcc = dash_module.dcc
    html = dash_module.html

    def _page_wrapper(page_id: str, content, visible: bool = False):
        style = {"display": "block"} if visible else {"display": "none"}
        return html.Div(content, id=page_id, style=style)

    return html.Div(
        [
            dcc.Location(id="url", refresh=False),
            dcc.Store(id="store-project", storage_type="local"),
            dcc.Store(id="store-run-config", storage_type="session"),
            dcc.Store(id="store-qc-outputs", storage_type="session"),
            dcc.Store(id="store-preprocess-outputs", storage_type="session"),
            dcc.Store(id="store-analysis-outputs", storage_type="session"),
            dcc.Store(id="store-databases-outputs", storage_type="session"),
            dcc.Store(id="store-phylo-outputs", storage_type="session"),
            dcc.Store(id="store-custom-db", storage_type="session"),
            dcc.Store(id="store-samples", storage_type="local"),
            html.Div(render_topbar(project_name="", mode="Complete", docker_ok=False), id="topbar", className="topbar"),
            html.Div(
                [
                    html.Div(render_sidebar(active_path="/project"), id="sidebar", className="sidebar"),
                    html.Div(
                        [
                            connection_banner(),
                            _page_wrapper("page-project", project_layout(), visible=True),
                            _page_wrapper("page-setup", setup_layout()),
                            _page_wrapper("page-samples", samples_layout()),
                            _page_wrapper("page-qc", qc_layout()),
                            _page_wrapper("page-preprocess", preprocess_layout()),
                            _page_wrapper("page-analysis", analysis_layout()),
                            _page_wrapper("page-databases", databases_layout()),
                            _page_wrapper("page-phylo", phylo_layout()),
                            _page_wrapper("page-progress", progress_layout()),
                            _page_wrapper("page-results", results_layout()),
                        ],
                        id="page-content",
                        className="main-content",
                    ),
                ],
                className="body-grid",
            ),
        ],
        className="app-shell",
    )


def get_layout():
    return create_layout()
