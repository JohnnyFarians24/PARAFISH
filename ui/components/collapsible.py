"""Collapsible output panels — for DETAIL, never for the headline.

Output is stratified into two layers, and the split is deliberate:

* **Layer 1, always visible** — the "did this run work?" layer: status badges, KPI
  cards, and ONE hero plot per page. A reviewer must be able to answer that question
  in three seconds, without clicking anything.
* **Layer 2, collapsed** — the audit layer: full tables, secondary plots, logs, raw
  headers, per-source dumps.

The rule is: *do I need this to trust the result?* → Layer 1. *Is this detail I go
looking for?* → Layer 2.

Collapsing everything (which an earlier pass did) is a real mistake in a scientific
tool: the page opens empty, the user cannot see whether the run succeeded, and hiding
the read-tracking and diversity tables actively undermines confidence in the numbers.

Collapsing is still preferred over making the user pick outputs *before* the run,
because at that point they cannot know what they will need — and re-running to recover
a chart they skipped is expensive (a remote BLAST is ~1 min per sequence). Only the
genuinely costly work is opt-in, on the Setup page.

``details()`` is the workhorse: native ``<details>``/``<summary>``, so the toggle costs
zero callbacks and zero JS. ``collapsible()`` uses ``dbc.Accordion`` and is for
several detail sections that should look like one stack.
"""

from __future__ import annotations

from typing import Optional

import dash_bootstrap_components as dbc
from dash import html


def details(title: str, children, *, subtitle: Optional[str] = None, open_by_default: bool = False):
    """A Layer-2 panel as a native ``<details>`` element. No callback, no JS."""
    summary = [html.Span(title, className="details-title")]
    if subtitle:
        summary.append(html.Span(subtitle, className="mini-note", style={"marginLeft": "10px"}))
    return html.Details(
        [
            html.Summary(summary, className="details-summary"),
            html.Div(children, className="details-body"),
        ],
        open=open_by_default,
        className="details-panel",
    )


def collapsible(title: str, children, *, subtitle: Optional[str] = None, item_id: Optional[str] = None, open_by_default: bool = False):
    """One collapsed output panel. Returns a ``dbc.Accordion`` holding a single item."""
    header = title if not subtitle else html.Span([html.Span(title), html.Span(subtitle, className="mini-note", style={"marginLeft": "10px", "opacity": 0.75})])
    item_id = item_id or title.lower().replace(" ", "-")
    return dbc.Accordion(
        [dbc.AccordionItem(children, title=header, item_id=item_id)],
        start_collapsed=not open_by_default,
        always_open=True,
        flush=True,
        active_item=[item_id] if open_by_default else [],
        className="output-accordion",
    )


def collapsible_group(panels: list[tuple[str, object]], *, open_first: bool = False):
    """Several output panels in one accordion, each independently openable.

    ``panels`` is a list of ``(title, content)``. ``always_open`` matters: opening one
    chart must not close the one the user was comparing it against.
    """
    items = []
    active: list[str] = []
    for index, (title, content) in enumerate(panels):
        item_id = f"panel-{index}-{str(title).lower().replace(' ', '-')}"
        items.append(dbc.AccordionItem(content, title=title, item_id=item_id))
        if index == 0 and open_first:
            active.append(item_id)
    return dbc.Accordion(items, start_collapsed=not open_first, always_open=True, flush=True, active_item=active, className="output-accordion")
