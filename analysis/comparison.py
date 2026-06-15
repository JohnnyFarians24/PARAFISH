"""Comparison structures for vsearch (OTU) vs dada2 (ASV) outputs."""

from __future__ import annotations

from datetime import datetime


def _normalize_summary(summary):
    return summary if isinstance(summary, dict) else {}


def _shared_keys(left, right):
    left_keys = set(left.keys())
    right_keys = set(right.keys())
    return sorted(left_keys.intersection(right_keys))


def _difference_keys(left, right):
    left_keys = set(left.keys())
    right_keys = set(right.keys())
    return {
        "vsearch_only": sorted(left_keys - right_keys),
        "dada2_only": sorted(right_keys - left_keys),
        "shared": _shared_keys(left, right),
    }


def compare_analysis_branches(vsearch_summary, dada2_summary):
    """Build a JSON-ready both-mode summary payload for analysis comparison."""
    vsearch = _normalize_summary(vsearch_summary)
    dada2 = _normalize_summary(dada2_summary)
    return {
        "mode": "both",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "vsearch": vsearch,
        "dada2": dada2,
        "comparison": _difference_keys(vsearch, dada2),
        "note": "ASV and OTU outputs are kept independent; comparison is derived from metrics only.",
    }
