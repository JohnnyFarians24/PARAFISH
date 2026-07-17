"""The single source of truth for how a DataTable looks.

Every page used to style its own tables. Three pages carried the same dark-theme dict
copy-pasted, and QC and Preprocess passed only ``{"padding": "8px"}`` — so those two
inherited Dash's default near-black cell colour and rendered dark-grey-on-dark: the
numbers were effectively invisible, and the only readable cells were the ones a
conditional rule happened to paint red. The page looked broken.

Styling tables in one place is the only way that divergence cannot come back. Call
``data_table()`` instead of ``dash_table.DataTable`` — there should be exactly one
``dash_table.DataTable(`` in the codebase, and it is in this file.
"""

from __future__ import annotations

from typing import Optional

from dash import dash_table

# #e5e7eb on the card background is ~13:1. The inherited default it replaces was ~2:1,
# well under the 4.5:1 floor.
TEXT = "#e5e7eb"
TEXT_MUTED = "#9ca3af"
RULE = "rgba(255,255,255,0.12)"

STYLE_TABLE = {"overflowX": "auto"}

STYLE_HEADER = {
    "backgroundColor": "transparent",
    "color": TEXT_MUTED,
    "fontWeight": 600,
    "border": "none",
    "borderBottom": f"1px solid {RULE}",
    "textAlign": "left",
}

STYLE_CELL = {
    "backgroundColor": "transparent",
    "color": TEXT,
    "border": "none",
    "fontFamily": "var(--font-mono, monospace)",
    "fontSize": "13px",
    "padding": "8px 12px",
    "textAlign": "left",
}

# Zebra striping does the work heavy cell borders used to, and keeps the eye on its row
# across a wide table.
STYLE_STRIPES = [{"if": {"row_index": "odd"}, "backgroundColor": "rgba(255,255,255,0.02)"}]

# Reliability / status wording is derived from identity + coverage, so colouring it makes
# the evidence legible at a glance instead of something to recompute per row.
SEMANTIC_STYLES = [
    {"if": {"filter_query": '{Reliability} = "high"', "column_id": "Reliability"}, "color": "#3fb950", "fontWeight": "600"},
    {"if": {"filter_query": '{Reliability} = "moderate"', "column_id": "Reliability"}, "color": "#d29922"},
    {"if": {"filter_query": '{Reliability} = "low"', "column_id": "Reliability"}, "color": "#f85149"},
    {"if": {"filter_query": '{reliability} = "high"', "column_id": "reliability"}, "color": "#3fb950", "fontWeight": "600"},
    {"if": {"filter_query": '{reliability} = "moderate"', "column_id": "reliability"}, "color": "#d29922"},
    {"if": {"filter_query": '{reliability} = "low"', "column_id": "reliability"}, "color": "#f85149"},
]

# Back-compat for callers that spread a style dict into DataTable(**TABLE_STYLE).
TABLE_STYLE = dict(style_table=STYLE_TABLE, style_header=STYLE_HEADER, style_cell=STYLE_CELL)


def data_table(data, columns, *, conditional: Optional[list] = None, semantic: bool = True, **kwargs):
    """A themed DataTable.

    ``conditional`` is *appended* to the base rules, never replaces them — a caller's
    red read-loss cells and the zebra striping have to coexist.
    """
    rules = [*STYLE_STRIPES, *(SEMANTIC_STYLES if semantic else []), *(conditional or [])]
    kwargs.setdefault("style_as_list_view", True)
    kwargs.setdefault("style_table", STYLE_TABLE)
    kwargs.setdefault("style_header", STYLE_HEADER)
    kwargs.setdefault("style_cell", STYLE_CELL)
    return dash_table.DataTable(
        data=data,
        columns=columns,
        style_data_conditional=rules,
        **kwargs,
    )
