"""Plot payload builders for analysis outputs."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go


def _build_heatmap_abundance(otu_matrix):
    if otu_matrix is None or otu_matrix.empty:
        return None

    normalized = otu_matrix.div(otu_matrix.sum(axis=0).replace(0, np.nan), axis=1).fillna(0) * 1_000_000.0
    return np.log1p(normalized)


def _build_stacked_abundance_table(otu_matrix, top_n=20):
    if otu_matrix is None or otu_matrix.empty:
        return None, None

    matrix = otu_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
    totals = matrix.sum(axis=1).sort_values(ascending=False)
    keep_features = list(totals.head(max(1, int(top_n))).index)
    selected = matrix.loc[keep_features].copy()
    other = matrix.drop(index=keep_features, errors="ignore")
    if not other.empty:
        selected.loc["Other"] = other.sum(axis=0)

    if selected.empty:
        return None, None

    long_rows = []
    for sample in selected.columns:
        sample_total = float(selected[sample].sum())
        for feature in selected.index:
            abundance = float(selected.at[feature, sample])
            long_rows.append(
                {
                    "sample": sample,
                    "feature": feature,
                    "abundance": abundance,
                    "relative_abundance": (abundance / sample_total) if sample_total else 0.0,
                    "sample_total": sample_total,
                }
            )

    long_df = pd.DataFrame(long_rows)
    fig = go.Figure()
    for feature in selected.index:
        fig.add_trace(
            go.Bar(
                name=str(feature),
                x=list(selected.columns),
                y=selected.loc[feature].astype(float).tolist(),
            )
        )
    fig.update_layout(
        barmode="stack",
        title="Stacked abundance by sample",
        xaxis_title="Sample",
        yaxis_title="Reads",
        legend_title="Feature",
        template="plotly_white",
    )
    return long_df, fig


def _build_rarefaction_curves(otu_matrix, max_points=12, replicates=5):
    if otu_matrix is None or otu_matrix.empty:
        return None, None

    matrix = otu_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
    rows = []
    fig = go.Figure()

    for sample in matrix.columns:
        counts = matrix[sample].astype(int).clip(lower=0)
        feature_labels = []
        for feature, count in counts.items():
            if count > 0:
                feature_labels.extend([str(feature)] * int(count))
        total_reads = len(feature_labels)
        if total_reads == 0:
            continue

        depths = np.unique(np.linspace(1, total_reads, num=min(int(max_points), total_reads), dtype=int))
        curve_replicates = []
        seed = abs(hash(sample)) % (2**32)
        rng = np.random.default_rng(seed)

        for _ in range(max(1, int(replicates))):
            permuted = np.array(feature_labels, dtype=object)
            rng.shuffle(permuted)
            seen = set()
            richness_values = []
            cursor = 0
            for depth in depths:
                while cursor < depth and cursor < total_reads:
                    seen.add(permuted[cursor])
                    cursor += 1
                richness_values.append(float(len(seen)))
            curve_replicates.append(richness_values)

        curve_array = np.array(curve_replicates, dtype=float)
        mean_curve = curve_array.mean(axis=0).tolist()
        std_curve = curve_array.std(axis=0).tolist() if len(curve_replicates) > 1 else [0.0] * len(mean_curve)

        for depth, richness_mean, richness_std in zip(depths, mean_curve, std_curve):
            rows.append(
                {
                    "sample": sample,
                    "depth": int(depth),
                    "richness_mean": round(float(richness_mean), 6),
                    "richness_std": round(float(richness_std), 6),
                    "total_reads": int(total_reads),
                }
            )

        fig.add_trace(
            go.Scatter(
                x=depths.tolist(),
                y=mean_curve,
                mode="lines+markers",
                name=str(sample),
            )
        )

    if not rows:
        return None, None

    fig.update_layout(
        title="Rarefaction curves",
        xaxis_title="Sequencing depth",
        yaxis_title="Observed richness",
        template="plotly_white",
    )
    return pd.DataFrame(rows), fig


def build_analysis_plot_data(comparison_data):
    """Build lightweight plotting payloads from analysis comparison data."""
    return comparison_data or {}
