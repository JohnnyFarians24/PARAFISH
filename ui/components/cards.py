"""Shared card, KPI and section-header components.

Phylo and Results drifted: they hand-rolled their own headings and stat blocks, so they
read as a different product from QC/Preprocess/Analysis. Everything visual that appears
on more than one page belongs here.
"""

from __future__ import annotations

from typing import Optional

from dash import html

# Status vocabulary used across every stage. One palette, one set of words.
STATUS_COLORS = {
    "success": ("#3fb950", "rgba(63,185,80,0.12)"),
    "partial": ("#d29922", "rgba(210,153,34,0.12)"),
    "failure": ("#f85149", "rgba(248,81,73,0.12)"),
    "skipped": ("#8b98aa", "rgba(139,152,170,0.12)"),
    "pending": ("#8b98aa", "rgba(139,152,170,0.12)"),
    "running": ("#58a6ff", "rgba(88,166,255,0.12)"),
}

STATUS_LABELS = {
    "success": "complete",
    "partial": "partial",
    "failure": "failed",
    "skipped": "not run",
    "pending": "pending",
    "running": "running",
}


def status_badge(status: str, *, label: Optional[str] = None):
    """A coloured pill naming a stage's outcome. Never renders blank."""
    key = str(status or "pending").lower()
    color, background = STATUS_COLORS.get(key, STATUS_COLORS["pending"])
    text = label or STATUS_LABELS.get(key, key)
    return html.Span(
        text,
        style={
            "color": color,
            "background": background,
            "border": f"1px solid {color}55",
            "borderRadius": "999px",
            "padding": "2px 10px",
            "fontSize": "11px",
            "fontWeight": "600",
            "whiteSpace": "nowrap",
        },
    )


def section(title: str, children, *, subtitle: Optional[str] = None, **kwargs):
    """A card with the accent-tick heading every page uses."""
    header = [html.Div(title, className="card-title")]
    if subtitle:
        header.append(html.Div(subtitle, className="mini-note"))
    return html.Div([*header, children], className="card card-body", **kwargs)


def kpi(value, label: str, *, hint: Optional[str] = None, color: Optional[str] = None, href: Optional[str] = None):
    """One KPI card: the single number that matters, plus what it means.

    ``href`` makes the whole card a link to the page that owns the number, so the
    Results summary can hand the reader straight to the detail.
    """
    body = [
        html.Div(str(value), className="stat-value", style={"color": color} if color else {}),
        html.Div(label, className="stat-label"),
    ]
    if hint:
        body.append(html.Div(hint, className="mini-note"))
    card = html.Div(body, className="stat-card")
    if not href:
        return card
    return html.A(card, href=href, style={"textDecoration": "none", "color": "inherit", "display": "block"})


def empty_state(message: str, *, hint: Optional[str] = None):
    """An explicit "nothing here, and here is why".

    A stage that did not run must say so. A blank panel is indistinguishable from a
    broken one, and in a scientific tool that difference matters.
    """
    body = [html.Div(message, className="empty-state-sub")]
    if hint:
        body.append(html.Div(hint, className="mini-note"))
    return html.Div(body, style={"padding": "12px 0"})
