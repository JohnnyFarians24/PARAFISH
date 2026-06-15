"""Centralized database selection logic by marker and source."""

from __future__ import annotations

from typing import Optional

from .bold_db import get_bold_reference
from .custom_db import get_custom_reference
from .ncbi_db import get_ncbi_reference
from .pr2_db import get_pr2_reference
from .silva_db import get_silva_reference


DATABASE_BY_SOURCE = {
    "ncbi": get_ncbi_reference,
    "silva": get_silva_reference,
    "pr2": get_pr2_reference,
    "bold": get_bold_reference,
    "custom": get_custom_reference,
}

DEFAULT_DATABASES_BY_MARKER = {
    "COI": ("bold", "ncbi", "custom"),
    "18S": ("silva", "pr2", "ncbi", "custom"),
    "ITS": ("silva", "pr2", "custom"),
    "28S": ("silva", "pr2", "custom"),
}


def resolve_database_references(marker: str, source: Optional[str] = None) -> list[str]:
    """Resolve one or many database references for a marker."""
    marker_key = marker.upper()
    if source:
        resolver = DATABASE_BY_SOURCE[source.lower()]
        return [resolver(marker_key)]

    allowed_sources = DEFAULT_DATABASES_BY_MARKER.get(marker_key, ("custom",))
    return [DATABASE_BY_SOURCE[name](marker_key) for name in allowed_sources]
