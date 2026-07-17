"""Cross-sample comparison metrics from an OTU/ASV abundance matrix.

The pipeline produces a feature x sample abundance table for both branches:
VSEARCH writes ``otu_matrix.csv`` (OTU rows, sample columns) and DADA2 writes
``asv_abundance.csv`` (ASV rows, sample columns). This module turns either into
the numbers a user needs to compare samples: per-sample richness and Shannon
diversity, how many features are shared across all samples versus exclusive to
one, and a compact heatmap payload.

It deliberately avoids pandas/plotly so it stays pure and unit-testable; the UI
layer turns the returned dict into figures.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Optional


def _load_matrix(path: Path) -> tuple[list[str], list[str], list[list[float]]]:
    """Load a feature x sample table.

    The first column is treated as the feature id; every remaining column whose
    values are numeric is treated as a sample. Returns
    ``(feature_ids, samples, values)`` where ``values[i][j]`` is the abundance of
    feature ``i`` in sample ``j``.
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        return [], [], []
    # Sniff delimiter (VSEARCH otu tables can be tab-separated; our CSVs comma).
    sample = text[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        delimiter = dialect.delimiter
    except Exception:
        delimiter = "\t" if "\t" in sample.splitlines()[0] else ","
    rows = list(csv.reader(text.splitlines(), delimiter=delimiter))
    rows = [r for r in rows if r]
    if len(rows) < 2:
        return [], [], []

    header = rows[0]
    body = rows[1:]

    # Decide which columns are numeric sample columns (skip the first / id column).
    candidate_cols = list(range(1, len(header)))

    def _num(value: str) -> Optional[float]:
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    numeric_cols = []
    for col in candidate_cols:
        seen = False
        ok = True
        for r in body:
            if col >= len(r):
                continue
            cell = r[col].strip()
            if cell == "":
                continue
            if _num(cell) is None:
                ok = False
                break
            seen = True
        if ok and seen:
            numeric_cols.append(col)

    samples = [header[col].strip() or f"sample_{col}" for col in numeric_cols]
    feature_ids: list[str] = []
    values: list[list[float]] = []
    for r in body:
        feature_ids.append((r[0].strip() if r else "") or f"feature_{len(feature_ids) + 1}")
        row_vals = []
        for col in numeric_cols:
            cell = r[col].strip() if col < len(r) else ""
            parsed = _num(cell)
            row_vals.append(parsed if parsed is not None else 0.0)
        values.append(row_vals)
    return feature_ids, samples, values


def _shannon(counts: list[float]) -> float:
    total = sum(c for c in counts if c > 0)
    if total <= 0:
        return 0.0
    h = 0.0
    for c in counts:
        if c > 0:
            p = c / total
            h -= p * math.log(p)
    return h


def build_comparison(path, top_n: int = 30) -> dict:
    """Build a cross-sample comparison payload from an abundance matrix file.

    Returns a dict with ``samples``, per-sample ``diversity`` (richness, Shannon,
    reads), ``shared``/``exclusive`` feature counts and a ``heatmap`` payload
    (top ``top_n`` features by total abundance). On any problem returns
    ``{"error": ...}`` with empty collections so the UI can degrade gracefully.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {"error": "Abundance matrix not found or empty", "samples": [], "diversity": [], "heatmap": None}

    try:
        feature_ids, samples, values = _load_matrix(path)
    except Exception as exc:  # pragma: no cover - defensive
        return {"error": f"Could not parse abundance matrix: {exc}", "samples": [], "diversity": [], "heatmap": None}

    if not samples or not feature_ids:
        return {"error": "No sample columns detected in abundance matrix", "samples": [], "diversity": [], "heatmap": None}

    n_samples = len(samples)

    # Per-sample diversity.
    diversity = []
    for j, sample in enumerate(samples):
        counts = [values[i][j] for i in range(len(feature_ids))]
        richness = sum(1 for c in counts if c > 0)
        diversity.append(
            {
                "sample": sample,
                "richness": richness,
                "shannon": round(_shannon(counts), 4),
                "reads": int(sum(counts)),
            }
        )

    # Shared vs exclusive: how many samples each feature is present in.
    shared = 0
    exclusive = {sample: 0 for sample in samples}
    for i in range(len(feature_ids)):
        present = [j for j in range(n_samples) if values[i][j] > 0]
        if len(present) == n_samples and n_samples > 0:
            shared += 1
        elif len(present) == 1:
            exclusive[samples[present[0]]] += 1

    # Heatmap: top features by total abundance across samples.
    totals = [(sum(values[i]), i) for i in range(len(feature_ids))]
    totals.sort(reverse=True)
    top = [i for _, i in totals[: max(1, top_n)]]
    heatmap = {
        "features": [feature_ids[i][:24] for i in top],
        "samples": samples,
        "values": [values[i] for i in top],
    }

    return {
        "samples": samples,
        "n_features": len(feature_ids),
        "n_samples": n_samples,
        "diversity": diversity,
        "shared": shared,
        "exclusive": exclusive,
        "heatmap": heatmap,
    }
