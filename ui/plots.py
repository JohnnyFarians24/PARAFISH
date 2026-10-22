"""Theme-aware Plotly figure builders for in-app output visualisation.

Figures are styled to sit on the PARAFISH dark cards: transparent backgrounds,
light fonts and the palette from ``assets/parafish_theme.css``. Each builder
returns a ``plotly.graph_objects.Figure`` (or ``None`` when there is nothing to
plot) so callbacks can wrap them in ``dcc.Graph`` and skip empties.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from analysis.ordination import DISTANCES
from phylo.phylo_metrics import SUPPORT_THRESHOLD
from preprocess.preprocess_metrics import MIN_STEP_RETENTION

if TYPE_CHECKING:  # pragma: no cover - typing only
    import pandas as pd

# Single source of truth for the "a step must not lose most of the reads" rule.
MIN_RETENTION_PCT = MIN_STEP_RETENTION * 100


# PARAFISH palette (mirrors assets/parafish_theme.css tokens).
_GREEN = "#3fb950"
_BLUE = "#67a6ff"
_TEAL = "#39d6b5"
_AMBER = "#d29922"
_PURPLE = "#a371f7"
_RED = "#f85149"
_SAND = "#d2a679"
_TEXT = "#c2ccd8"
_GRID = "rgba(139,152,170,0.15)"

# Categorical palette for taxa (enough distinct hues for top-N + Other).
_CATEGORICAL = [
    _BLUE, _TEAL, _GREEN, _AMBER, _PURPLE, _SAND, _RED,
    "#8bb9ff", "#5ad1b6", "#7ee08a", "#e8c56b", "#c39bf5", "#e0b892", "#ff8a82",
    "#8b98aa",  # grey reserved for the "Other" bucket tail
]


def _theme(fig: go.Figure, height: int = 320, title: Optional[str] = None) -> go.Figure:
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color=_TEXT, size=12),
        title=dict(text=title or "", font=dict(color="#f4f7fb", size=14)),
        margin=dict(l=48, r=16, t=40 if title else 16, b=44),
        height=height,
        legend=dict(bgcolor="rgba(0,0,0,0)", orientation="h", yanchor="bottom", y=1.02, x=0),
        colorway=[_BLUE, _TEAL, _GREEN, _AMBER, _PURPLE],
        # Space between bars, so grouped/stacked charts read as bars and not as one solid
        # wall. Ignored by non-bar figures, so this is safe to apply to every figure.
        bargap=0.32,
        bargroupgap=0.12,
    )
    fig.update_xaxes(gridcolor=_GRID, zerolinecolor=_GRID)
    fig.update_yaxes(gridcolor=_GRID, zerolinecolor=_GRID)
    return fig


# ── Cross-sample comparison ─────────────────────────────────────────────

def comparison_figures(data: dict) -> list[go.Figure]:
    """Build the cross-sample comparison figures from ``build_comparison`` data."""
    figures: list[go.Figure] = []
    if not data or data.get("error") or not data.get("samples"):
        return figures

    # 1. Abundance heatmap (top features x samples).
    heatmap = data.get("heatmap")
    if heatmap and heatmap.get("features") and heatmap.get("samples"):
        fig = go.Figure(
            go.Heatmap(
                z=heatmap["values"],
                x=heatmap["samples"],
                y=heatmap["features"],
                colorscale="Viridis",
                colorbar=dict(title="reads"),
                hovertemplate="Feature %{y}<br>Sample %{x}<br>%{z} reads<extra></extra>",
            )
        )
        fig.update_layout(yaxis=dict(autorange="reversed"))
        figures.append(_theme(fig, height=max(300, 18 * len(heatmap["features"]) + 120), title="Abundance heatmap (top features × sample)"))

    # 2. Per-sample diversity: richness (bars) + Shannon (line, secondary axis).
    diversity = data.get("diversity") or []
    if diversity:
        samples = [d["sample"] for d in diversity]
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Bar(x=samples, y=[d["richness"] for d in diversity], name="Richness (features)", marker_color=_BLUE), secondary_y=False)
        fig.add_trace(go.Scatter(x=samples, y=[d["shannon"] for d in diversity], name="Shannon", mode="lines+markers", marker_color=_AMBER), secondary_y=True)
        _theme(fig, height=320, title="Per-sample diversity")
        fig.update_yaxes(title_text="Richness", secondary_y=False, gridcolor=_GRID)
        fig.update_yaxes(title_text="Shannon", secondary_y=True, showgrid=False)
        figures.append(fig)

    # 3. Shared vs exclusive features.
    exclusive = data.get("exclusive") or {}
    if exclusive:
        samples = list(exclusive.keys())
        fig = go.Figure()
        fig.add_trace(go.Bar(x=samples, y=[exclusive[s] for s in samples], name="Exclusive to sample", marker_color=_TEAL))
        shared = data.get("shared", 0)
        fig.add_hline(y=shared, line_dash="dash", line_color=_GREEN, annotation_text=f"Shared by all: {shared}", annotation_font_color=_GREEN)
        figures.append(_theme(fig, height=300, title="Exclusive features per sample (dashed = shared by all samples)"))

    return figures


# ── QC ──────────────────────────────────────────────────────────────────
#
# Total reads and %GC are deliberately *not* plotted: both already appear in the
# Run overview table, and a bar chart of a column you can already read is noise.
# What the table cannot show is a distribution or a per-position profile.

# Below this many samples a distribution plot is a lie: a box drawn through two points
# shows a quartile range that does not exist, and the eye reads the empty space as
# missing data. A labelled bar per sample says exactly what is known and nothing more.
SMALL_N = 4


def qc_depth_box_figure(overview: "pd.DataFrame") -> Optional[go.Figure]:
    """Read depth across samples: a bar per sample when few, a distribution when many.

    Uses R1 only: in PE data R2 carries the same read count, so including it would
    double every point and halve the apparent spread.
    """
    r1 = overview[overview["read"] == "R1"].dropna(subset=["total_reads"])
    if r1.empty:
        return None
    outliers = r1[r1["is_depth_outlier"] == 1]

    if len(r1) <= SMALL_N:
        colors = [_RED if flag == 1 else _BLUE for flag in r1["is_depth_outlier"]]
        fig = go.Figure(
            go.Bar(
                x=r1["sample"],
                y=r1["total_reads"],
                marker_color=colors,
                text=[f"{value:,.0f}" for value in r1["total_reads"]],
                textposition="outside",
                hovertemplate="%{x}<br>%{y:,.0f} reads<extra></extra>",
            )
        )
        title = f"Read depth — {len(r1)} sample(s)"
        if not outliers.empty:
            title += f" — {len(outliers)} low-depth outlier(s) in red"
        fig = _theme(fig, height=320, title=title)
        fig.update_yaxes(title_text="Reads (R1)")
        return fig

    if len(r1) < 2:
        return None
    fig = go.Figure(
        go.Box(
            y=r1["total_reads"],
            name="R1 reads",
            boxpoints="all",
            jitter=0.4,
            marker_color=_BLUE,
            line_color=_TEAL,
            text=r1["sample"],
            hovertemplate="%{text}<br>%{y:,.0f} reads<extra></extra>",
        )
    )
    if not outliers.empty:
        fig.add_trace(
            go.Scatter(
                x=[0] * len(outliers),
                y=outliers["total_reads"],
                mode="markers",
                name="Low-depth outlier",
                marker=dict(color=_RED, size=11, symbol="x"),
                text=outliers["sample"],
                hovertemplate="%{text}<br>%{y:,.0f} reads — below half the median<extra></extra>",
            )
        )
    title = f"Read-depth distribution across {len(r1)} samples"
    if not outliers.empty:
        title += f" — {len(outliers)} low-depth outlier(s)"
    return _theme(fig, height=320, title=title)


def qc_per_base_quality_figure(per_base: "pd.DataFrame") -> Optional[go.Figure]:
    """Mean Phred score per position, one line per sample x read."""
    if per_base.empty:
        return None
    fig = go.Figure()
    for (sample, read), group in per_base.groupby(["sample", "read"], sort=True):
        group = group.sort_values("position")
        fig.add_trace(
            go.Scatter(
                x=group["position"],
                y=group["mean_quality"],
                mode="lines",
                name=f"{sample} {read}",
                line=dict(width=1.5, dash="solid" if read == "R1" else "dot"),
            )
        )
    # Q20 = 1 error in 100 bases, Q30 = 1 in 1000: the conventional accept/reject bands.
    fig.add_hrect(y0=0, y1=20, fillcolor=_RED, opacity=0.08, line_width=0)
    fig.add_hrect(y0=20, y1=30, fillcolor=_AMBER, opacity=0.08, line_width=0)
    fig.add_hline(y=30, line_dash="dash", line_color=_GREEN, annotation_text="Q30", annotation_font_color=_GREEN)
    _theme(fig, height=360, title="Mean quality per base position")
    fig.update_xaxes(title_text="Position in read (bp)")
    fig.update_yaxes(title_text="Mean Phred score", rangemode="tozero")
    return fig


def qc_length_distribution_figure(lengths: "pd.DataFrame") -> Optional[go.Figure]:
    """Read-length distribution per sample, as overlaid density-style lines."""
    if lengths.empty:
        return None
    # A single fixed length (untrimmed Illumina output) plots as one spike per
    # sample, which reads better as bars than as a line going nowhere.
    single_bin = lengths.groupby(["sample", "read"])["length"].nunique().max() == 1
    fig = go.Figure()
    for (sample, read), group in lengths.groupby(["sample", "read"], sort=True):
        group = group.sort_values("length")
        name = f"{sample} {read}"
        if single_bin:
            fig.add_trace(go.Bar(x=group["length"], y=group["count"], name=name))
        else:
            fig.add_trace(go.Scatter(x=group["length"], y=group["count"], mode="lines", name=name, line=dict(width=1.5)))
    _theme(fig, height=340, title="Read-length distribution")
    fig.update_xaxes(title_text="Read length (bp)")
    fig.update_yaxes(title_text="Reads")
    return fig


# ── Preprocess ──────────────────────────────────────────────────────────

def read_fate_figure(tracking: "pd.DataFrame") -> Optional[go.Figure]:
    """Where every input read ended up, stacked per sample.

    The categories partition the input exactly, so a bar that does not reach the
    input count means the accounting is broken, not that reads vanished quietly.
    """
    import pandas as pd

    if tracking.empty:
        return None
    samples = tracking["sample"].tolist()
    categories = [
        ("assembled_reads", "Merged (kept)", _GREEN),
        ("discarded_reads", "Discarded by PEAR", _AMBER),
        ("not_assembled_reads", "Not merged", _PURPLE),
        ("forward_only", "Forward only", _BLUE),
        ("reverse_only", "Reverse only", _TEAL),
        ("dropped_reads", "Dropped by trimming", _RED),
    ]
    fig = go.Figure()
    for column, label, colour in categories:
        if column not in tracking:
            continue
        # Single-end rows leave the PEAR columns as empty strings, and "" is neither
        # NaN nor 0 — without coercing, those categories are drawn as empty traces.
        values = pd.to_numeric(tracking[column], errors="coerce").fillna(0)
        if values.sum() == 0:
            continue
        fig.add_trace(go.Bar(x=samples, y=values, name=label, marker_color=colour))
    if not fig.data:
        return None
    fig.update_layout(barmode="stack")
    _theme(fig, height=340, title="Read fate per sample")
    fig.update_yaxes(title_text="Read pairs")
    return fig


def assembly_rate_figure(tracking: "pd.DataFrame") -> Optional[go.Figure]:
    """Share of trimmed pairs PEAR managed to merge, per sample.

    Returns ``None`` for single-end runs, where there is nothing to merge.
    """
    import pandas as pd

    if "pct_assembled" not in tracking:
        return None
    values = pd.to_numeric(tracking["pct_assembled"], errors="coerce")
    keep = values.notna()
    if not keep.any():
        return None
    rates = tracking[keep]
    values = values[keep]
    colours = [_RED if value < MIN_RETENTION_PCT else _GREEN for value in values]
    fig = go.Figure(
        go.Bar(
            x=rates["sample"],
            y=values,
            marker_color=colours,
            hovertemplate="%{x}<br>%{y:.2f}% of trimmed pairs merged<extra></extra>",
        )
    )
    fig.add_hline(
        y=MIN_RETENTION_PCT,
        line_dash="dash",
        line_color=_RED,
        annotation_text=f"{MIN_RETENTION_PCT:.0f}% — below this, most pairs fail to merge",
        annotation_font_color=_RED,
    )
    _theme(fig, height=300, title="Assembly rate (paired-end merging)")
    fig.update_yaxes(title_text="% merged", range=[0, 100])
    return fig


def fragment_length_figure(lengths: "pd.DataFrame") -> Optional[go.Figure]:
    """Merged-fragment length distribution: a clean peak means the amplicon merged well."""
    if lengths.empty:
        return None
    fig = go.Figure()
    for sample, group in lengths.groupby("sample", sort=True):
        group = group.sort_values("length")
        fig.add_trace(go.Scatter(x=group["length"], y=group["count"], mode="lines", name=str(sample), line=dict(width=1.5)))
    _theme(fig, height=340, title="Merged fragment length distribution")
    fig.update_xaxes(title_text="Fragment length (bp)")
    fig.update_yaxes(title_text="Reads")
    return fig


# ── Analysis: diversity ─────────────────────────────────────────────────

def rarefaction_figure(curve: "pd.DataFrame", depth: Optional[int] = None) -> Optional[go.Figure]:
    """Richness against sequencing depth. Curves that flatten mean enough was sequenced."""
    if curve.empty:
        return None
    fig = go.Figure()
    for sample, group in curve.groupby("sample", sort=True):
        group = group.sort_values("depth")
        fig.add_trace(go.Scatter(x=group["depth"], y=group["richness"], mode="lines", name=str(sample), line=dict(width=1.6)))
    if depth:
        fig.add_vline(
            x=depth,
            line_dash="dash",
            line_color=_AMBER,
            annotation_text=f"rarefaction depth {depth:,}",
            annotation_font_color=_AMBER,
        )
    _theme(fig, height=360, title="Rarefaction curves")
    fig.update_xaxes(title_text="Reads sampled")
    fig.update_yaxes(title_text="Observed features")
    return fig


def rank_abundance_figure(ranks: "pd.DataFrame") -> Optional[go.Figure]:
    """Whittaker plot: a steep curve means a few features dominate the sample."""
    if ranks.empty:
        return None
    fig = go.Figure()
    for sample, group in ranks.groupby("sample", sort=True):
        group = group.sort_values("rank")
        fig.add_trace(
            go.Scatter(x=group["rank"], y=group["relative_abundance"], mode="lines+markers", name=str(sample), marker=dict(size=4))
        )
    _theme(fig, height=340, title="Rank-abundance (Whittaker)")
    fig.update_xaxes(title_text="Rank")
    fig.update_yaxes(title_text="Relative abundance", type="log")
    return fig


def alpha_diversity_figure(alpha: "pd.DataFrame") -> Optional[go.Figure]:
    """Alpha diversity across samples, one panel per metric.

    With few samples this draws a labelled bar per sample instead of a box: a box
    through two points invents a spread, and the reader wants the per-sample value
    anyway. With many samples the distribution is the point, so the box returns.
    """
    if alpha.empty:
        return None
    metrics = [("observed", "Observed richness"), ("shannon", "Shannon"), ("pielou", "Pielou evenness")]
    small = len(alpha) <= SMALL_N
    fig = make_subplots(rows=1, cols=len(metrics), subplot_titles=[label for _, label in metrics])
    for index, (column, _) in enumerate(metrics, start=1):
        if column not in alpha:
            continue
        if small:
            trace = go.Bar(
                x=alpha["sample"],
                y=alpha[column],
                marker_color=_BLUE,
                text=[f"{value:.2f}" if isinstance(value, float) else str(value) for value in alpha[column]],
                textposition="outside",
                hovertemplate="%{x}<br>%{y}<extra></extra>",
                showlegend=False,
            )
        else:
            trace = go.Box(
                y=alpha[column],
                name="",
                boxpoints="all",
                jitter=0.5,
                pointpos=0,
                marker_color=_BLUE,
                line_color=_TEAL,
                text=alpha["sample"],
                hovertemplate="%{text}<br>%{y}<extra></extra>",
                showlegend=False,
            )
        fig.add_trace(trace, row=1, col=index)
    suffix = f" — {len(alpha)} sample(s)" if small else ""
    _theme(fig, height=320, title=f"Alpha diversity (rarefied){suffix}")
    return fig


def feature_length_figure(lengths: "pd.DataFrame") -> Optional[go.Figure]:
    """OTU/ASV length distribution. Peaks away from the amplicon size are artefacts."""
    if lengths.empty:
        return None
    fig = go.Figure(go.Bar(x=lengths["length"], y=lengths["count"], marker_color=_PURPLE))
    _theme(fig, height=300, title="OTU / ASV length distribution")
    fig.update_xaxes(title_text="Length (bp)")
    fig.update_yaxes(title_text="Features")
    return fig


def upset_figure(intersections: "pd.DataFrame", samples: list[str]) -> Optional[go.Figure]:
    """UpSet plot: intersection sizes above, a dot matrix of membership below.

    Replaces the Venn diagram, which stops being legible past three sets.
    """
    if intersections.empty or not samples:
        return None
    data = intersections.sort_values("n_features", ascending=False).head(20).reset_index(drop=True)
    positions = list(range(len(data)))
    memberships = [set(str(row).split("|")) for row in data["combination"]]

    fig = make_subplots(rows=2, cols=1, row_heights=[0.62, 0.38], shared_xaxes=True, vertical_spacing=0.06)
    fig.add_trace(
        go.Bar(
            x=positions,
            y=data["n_features"],
            marker_color=_BLUE,
            text=data["n_features"],
            textposition="outside",
            hovertemplate="%{customdata}<br>%{y} features<extra></extra>",
            customdata=data["combination"],
            showlegend=False,
        ),
        row=1,
        col=1,
    )
    # Dot matrix: a filled dot means the sample participates in that intersection.
    for row_index, sample in enumerate(samples):
        y = [row_index] * len(positions)
        filled = [index for index, members in enumerate(memberships) if sample in members]
        empty = [index for index in positions if index not in filled]
        if empty:
            fig.add_trace(
                go.Scatter(x=empty, y=[row_index] * len(empty), mode="markers", marker=dict(size=9, color="rgba(139,152,170,0.25)"), showlegend=False, hoverinfo="skip"),
                row=2,
                col=1,
            )
        if filled:
            fig.add_trace(
                go.Scatter(x=filled, y=[row_index] * len(filled), mode="markers", marker=dict(size=10, color=_TEAL), showlegend=False, name=sample, hovertemplate=f"{sample}<extra></extra>"),
                row=2,
                col=1,
            )
        del y
    _theme(fig, height=max(340, 90 + 26 * len(samples) + 180), title="Shared and exclusive features (UpSet)")
    fig.update_yaxes(title_text="Features", row=1, col=1)
    fig.update_yaxes(
        tickmode="array",
        tickvals=list(range(len(samples))),
        ticktext=samples,
        showgrid=False,
        row=2,
        col=1,
    )
    fig.update_xaxes(showticklabels=False, showgrid=False, row=1, col=1)
    fig.update_xaxes(showticklabels=False, showgrid=False, row=2, col=1)
    return fig


# ── Databases: taxonomic resolution and confidence ──────────────────────

_RANK_LABELS = {
    "kingdom": "Kingdom",
    "phylum": "Phylum",
    "class_": "Class",
    "order": "Order",
    "family": "Family",
    "genus": "Genus",
    "species": "Species",
    "unassigned": "Unassigned",
}


def taxonomic_resolution_figure(resolution: "pd.DataFrame") -> Optional[go.Figure]:
    """How deep the taxonomy got, counted by feature *and* by read.

    The two bars routinely disagree: many unresolved features can carry a sliver of
    the reads, or one unresolved feature can carry half the run.
    """
    if resolution.empty:
        return None
    labels = [_RANK_LABELS.get(rank, str(rank)) for rank in resolution["rank"]]
    fig = go.Figure()
    fig.add_trace(go.Bar(x=labels, y=resolution["pct_features"], name="% of features", marker_color=_BLUE))
    fig.add_trace(go.Bar(x=labels, y=resolution["pct_reads"], name="% of reads", marker_color=_TEAL))
    fig.update_layout(barmode="group")
    _theme(fig, height=340, title="Taxonomic resolution — deepest rank reached")
    fig.update_yaxes(title_text="%")
    return fig


def confidence_histogram_figure(assignments: "pd.DataFrame") -> Optional[go.Figure]:
    """Distributions of identity, coverage and bootstrap across all assignments.

    Panels are omitted when the sources in this run do not produce that evidence:
    SINTAX never aligns, so an identity panel over a SINTAX-only run would be an
    empty axis pretending a measurement exists.
    """
    import pandas as pd

    if assignments.empty:
        return None
    panels = []
    for column, label in (("identity", "% identity"), ("coverage", "% query coverage"), ("confidence", "Confidence / bootstrap")):
        values = pd.to_numeric(assignments.get(column), errors="coerce").dropna() if column in assignments else pd.Series(dtype=float)
        if not values.empty:
            panels.append((label, values))
    if not panels:
        return None

    fig = make_subplots(rows=1, cols=len(panels), subplot_titles=[label for label, _ in panels])
    for index, (_, values) in enumerate(panels, start=1):
        fig.add_trace(go.Histogram(x=values, marker_color=_PURPLE, showlegend=False), row=1, col=index)
    _theme(fig, height=300, title="Assignment confidence")
    return fig


def reliability_figure(assignments: "pd.DataFrame") -> Optional[go.Figure]:
    """Reliability of each source's calls, per rank."""
    if assignments.empty or "reliability" not in assignments:
        return None
    counts = assignments.groupby(["source", "reliability"]).size().reset_index(name="n")
    if counts.empty:
        return None
    colours = {"high": _GREEN, "moderate": _AMBER, "low": _RED, "unknown": "#8b98aa", "unassigned": "#5a6472"}
    fig = go.Figure()
    for reliability, group in counts.groupby("reliability"):
        fig.add_trace(
            go.Bar(x=group["source"], y=group["n"], name=str(reliability), marker_color=colours.get(str(reliability), _PURPLE))
        )
    fig.update_layout(barmode="stack")
    _theme(fig, height=300, title="Assignment reliability by source")
    fig.update_yaxes(title_text="Assignments")
    return fig


def taxonomy_sunburst_figure(consensus: "pd.DataFrame", max_depth: int = 6) -> Optional[go.Figure]:
    """Interactive hierarchical view of the consensus taxonomy.

    Stands in for Krona, which is a Perl toolchain rather than a pip package: a
    sunburst gives the same click-to-zoom exploration with no extra dependency.
    """
    columns = ["Consensus_Kingdom", "Consensus_Phylum", "Consensus_Class", "Consensus_Order", "Consensus_Family", "Consensus_Genus"]
    present = [column for column in columns if column in consensus]
    if consensus.empty or not present:
        return None

    frame = consensus[present].fillna("").astype(str)
    frame = frame.replace({"": "Unassigned", "nan": "Unassigned"})

    labels: list[str] = []
    parents: list[str] = []
    values: list[int] = []
    index_of: dict[str, int] = {}
    for _, row in frame.iterrows():
        path = ""
        for column in present[:max_depth]:
            parent = path
            path = f"{path}/{row[column]}" if path else row[column]
            if path not in index_of:
                index_of[path] = len(labels)
                labels.append(path)
                parents.append(parent)
                values.append(0)
            # Each ancestor accumulates one per descendant leaf, so a parent's value
            # is exactly the sum of its children — what branchvalues="total" requires.
            values[index_of[path]] += 1

    if not labels:
        return None
    display = [label.split("/")[-1] for label in labels]
    fig = go.Figure(
        go.Sunburst(ids=labels, labels=display, parents=parents, values=values, branchvalues="total", hovertemplate="%{id}<br>%{value} features<extra></extra>")
    )
    return _theme(fig, height=520, title="Taxonomic composition (click to zoom)")


def prevalence_figure(prevalence: "pd.DataFrame", top_n: int = 20) -> Optional[go.Figure]:
    """Per-taxon prevalence with 95% Wilson intervals.

    The interval is the point: a taxon seen in 3 of 5 hosts is 60% prevalent with
    a 95% interval from 23% to 88%, and the bar alone hides that.
    """
    if prevalence.empty:
        return None
    data = prevalence.head(top_n).iloc[::-1]
    fig = go.Figure(
        go.Bar(
            x=data["prevalence"] * 100,
            y=data["taxon"],
            orientation="h",
            marker_color=_TEAL,
            error_x=dict(
                type="data",
                symmetric=False,
                array=(data["ci_high"] - data["prevalence"]) * 100,
                arrayminus=(data["prevalence"] - data["ci_low"]) * 100,
                color=_TEXT,
            ),
            hovertemplate="%{y}<br>%{x:.1f}% (95% CI)<extra></extra>",
        )
    )
    _theme(fig, height=max(320, 22 * len(data) + 120), title="Taxon prevalence across samples (95% Wilson CI)")
    fig.update_xaxes(title_text="Prevalence (%)", range=[0, 100])
    return fig


def contaminant_figure(contaminants: "pd.DataFrame") -> Optional[go.Figure]:
    """Prevalence in the negative controls against prevalence in real samples."""
    if contaminants.empty:
        return None
    flagged = contaminants["is_contaminant"] == 1
    fig = go.Figure()
    for mask, name, colour in ((~flagged, "Kept", _TEAL), (flagged, "Flagged as contaminant", _RED)):
        subset = contaminants[mask]
        if subset.empty:
            continue
        fig.add_trace(
            go.Scatter(
                x=subset["prevalence_samples"] * 100,
                y=subset["prevalence_controls"] * 100,
                mode="markers",
                name=name,
                marker=dict(color=colour, size=9),
                text=subset["feature"],
                hovertemplate="%{text}<br>samples %{x:.0f}%<br>controls %{y:.0f}%<extra></extra>",
            )
        )
    if not fig.data:
        return None
    _theme(fig, height=380, title="Contamination: prevalence in controls vs samples")
    fig.update_xaxes(title_text="Prevalence in biological samples (%)", range=[-5, 105])
    fig.update_yaxes(title_text="Prevalence in negative controls (%)", range=[-5, 105])
    return fig


# ── Taxonomy / community ────────────────────────────────────────────────

def composition_figure(comp: dict) -> Optional[go.Figure]:
    """Stacked relative-abundance bar (samples on x, taxa as stacked segments)."""
    if not comp or comp.get("error") or not comp.get("samples") or not comp.get("taxa"):
        return None
    samples = comp["samples"]
    taxa = comp["taxa"]
    rel = comp["rel"]
    fig = go.Figure()
    for t, taxon in enumerate(taxa):
        color = "#8b98aa" if taxon == "Other" else _CATEGORICAL[t % (len(_CATEGORICAL) - 1)]
        fig.add_trace(
            go.Bar(
                x=samples,
                y=[rel[t][j] * 100 for j in range(len(samples))],
                name=taxon,
                marker_color=color,
                hovertemplate="%{x}<br>" + taxon + ": %{y:.1f}%<extra></extra>",
            )
        )
    fig.update_layout(barmode="stack", yaxis_title="Relative abundance (%)")
    height = 360 if len(samples) <= 40 else 420
    return _theme(fig, height=height, title=f"Taxonomic composition by {comp.get('rank', 'taxon')} (top taxa)")


def ordination_figure(pcoa_result: dict) -> Optional[go.Figure]:
    """PCoA ordination scatter — one point per sample; nearby = similar community."""
    if not pcoa_result or len(pcoa_result.get("samples", [])) < 3:
        return None
    samples = pcoa_result["samples"]
    coords = pcoa_result["coords"]
    var = pcoa_result.get("variance", [0, 0])
    fig = go.Figure(
        go.Scatter(
            x=[c[0] for c in coords],
            y=[c[1] for c in coords],
            mode="markers+text" if len(samples) <= 30 else "markers",
            text=samples if len(samples) <= 30 else None,
            textposition="top center",
            textfont=dict(size=9),
            marker=dict(size=11, color=_TEAL, line=dict(width=1, color="#0b0f14")),
            hovertext=samples,
            hovertemplate="%{hovertext}<br>PC1 %{x:.3f}<br>PC2 %{y:.3f}<extra></extra>",
        )
    )
    fig.update_layout(xaxis_title=f"PCoA 1 ({var[0]}%)", yaxis_title=f"PCoA 2 ({var[1]}%)")
    return _theme(fig, height=420, title="Sample ordination (Bray-Curtis PCoA)")


def beta_heatmap_figure(bc: dict) -> Optional[go.Figure]:
    """Sample x sample Bray-Curtis dissimilarity heatmap."""
    samples = (bc or {}).get("samples") or []
    matrix = (bc or {}).get("matrix") or []
    if len(samples) < 2 or len(samples) > 60:  # unreadable beyond ~60 labels
        return None
    fig = go.Figure(
        go.Heatmap(
            z=matrix, x=samples, y=samples, colorscale="Magma", zmin=0, zmax=1,
            colorbar=dict(title="Bray-Curtis"),
            hovertemplate="%{y} vs %{x}<br>dissimilarity %{z:.3f}<extra></extra>",
        )
    )
    fig.update_layout(yaxis=dict(autorange="reversed"))
    return _theme(fig, height=max(320, 16 * len(samples) + 120), title="Between-sample dissimilarity (Bray-Curtis)")


# ── Custom database comparison ──────────────────────────────────────────

def detection_heatmap_figure(detection: dict, *, relative: bool = True, max_taxa: int = 40) -> Optional[go.Figure]:
    """Taxon × sample detection heatmap for the custom-DB comparison.

    Columns are samples, rows are taxa (matched species + a ``No match`` row).
    Column-normalises to relative abundance by default so the colour scale stays
    comparable across dozens to hundreds of samples with very different depths.
    """
    if not detection:
        return None
    taxa = detection.get("taxa") or []
    samples = detection.get("samples") or []
    values = detection.get("values") or []
    if not taxa or not samples or not values:
        return None
    # Degenerate case: nothing matched, so the only row is "No match". A single
    # column-normalised row is 100% everywhere — a misleading uniform heatmap.
    # Skip it; the match-rate donut already communicates "0% matched".
    matched_taxa = [t for t in taxa if str(t).strip().lower() != "no match"]
    if not matched_taxa:
        return None

    # Keep the most abundant taxa readable; the caller already ordered them.
    if len(taxa) > max_taxa:
        taxa = taxa[:max_taxa]
        values = values[:max_taxa]

    z = [[float(v) for v in row] for row in values]
    unit = "reads"
    if relative:
        col_totals = [sum(z[t][j] for t in range(len(z))) for j in range(len(samples))]
        z = [[(z[t][j] / col_totals[j] * 100.0 if col_totals[j] > 0 else 0.0) for j in range(len(samples))] for t in range(len(z))]
        unit = "%"

    fig = go.Figure(
        go.Heatmap(
            z=z,
            x=samples,
            y=taxa,
            colorscale="Viridis",
            colorbar=dict(title=unit),
            hovertemplate="%{y}<br>%{x}<br>%{z:.1f} " + unit + "<extra></extra>",
        )
    )
    fig.update_layout(yaxis=dict(autorange="reversed"))
    height = max(300, 22 * len(taxa) + 130)
    return _theme(fig, height=height, title=f"Detection per sample ({'relative abundance' if relative else 'reads'}) — custom database comparison")


def match_rate_figure(summary: dict) -> Optional[go.Figure]:
    """Donut of matched vs unmatched OTUs against the custom reference."""
    if not summary:
        return None
    matched = int(summary.get("n_matched", 0))
    unmatched = int(summary.get("n_unmatched", 0))
    if matched + unmatched <= 0:
        return None
    fig = go.Figure(
        go.Pie(
            labels=["Matched", "None (no match)"],
            values=[matched, unmatched],
            hole=0.55,
            marker=dict(colors=[_GREEN, _RED]),
            textinfo="label+percent",
            hovertemplate="%{label}: %{value} OTUs (%{percent})<extra></extra>",
        )
    )
    return _theme(fig, height=300, title=f"OTUs matched in the custom database ({summary.get('pct_matched', 0)}%)")


# ── Phylogenetics ───────────────────────────────────────────────────────

def pd_figure(faith: dict) -> Optional[go.Figure]:
    """Per-sample Faith's phylogenetic diversity bar chart."""
    if not faith or faith.get("error") or not faith.get("samples"):
        return None
    fig = go.Figure(go.Bar(x=faith["samples"], y=faith["pd"], marker_color=_PURPLE, name="Faith's PD"))
    fig.update_layout(yaxis_title="Faith's PD (branch length)")
    return _theme(fig, height=300, title="Phylogenetic diversity per sample (Faith's PD)")


def read_fate_sankey_figure(sankey: dict) -> Optional[go.Figure]:
    """End-to-end read fate, from raw reads to classified reads, as a Sankey."""
    if not sankey or not sankey.get("stages"):
        return None
    stages = sankey["stages"]
    index = {name: position for position, name in enumerate(stages)}
    colours = [_RED, _AMBER, _PURPLE, _BLUE, _TEAL, _GREEN]
    node_labels = [f"{name}<br>{value:,}" for name, value in zip(stages, sankey["values"])]
    fig = go.Figure(
        go.Sankey(
            node=dict(
                label=node_labels,
                color=[colours[position % len(colours)] for position in range(len(stages))],
                pad=18,
                thickness=18,
                line=dict(width=0),
            ),
            link=dict(
                source=[index[flow["source"]] for flow in sankey["flows"]],
                target=[index[flow["target"]] for flow in sankey["flows"]],
                value=[flow["value"] for flow in sankey["flows"]],
                color="rgba(103,166,255,0.25)",
            ),
        )
    )
    return _theme(fig, height=360, title="Read fate across the pipeline")


def dashboard_richness_figure(rows: list[dict]) -> Optional[go.Figure]:
    """Reads and richness per sample, controls drawn apart from real samples."""
    if not rows:
        return None
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    samples = [row["sample"] for row in rows]
    colours = [_RED if row.get("is_control") else _BLUE for row in rows]
    fig.add_trace(
        go.Bar(x=samples, y=[row["reads"] for row in rows], name="Surviving reads", marker_color=colours),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(x=samples, y=[row["richness"] for row in rows], name="Richness", mode="lines+markers", marker_color=_AMBER),
        secondary_y=True,
    )
    _theme(fig, height=320, title="Per-sample reads and richness (controls in red)")
    fig.update_yaxes(title_text="Reads", secondary_y=False, gridcolor=_GRID)
    fig.update_yaxes(title_text="Richness", secondary_y=True, showgrid=False)
    return fig


def alignment_profile_figure(columns: "pd.DataFrame") -> Optional[go.Figure]:
    """Per-column conservation and gap fraction of the MAFFT alignment.

    Both curves are needed: a column that is 92% gaps but whose single remaining
    base is invariant scores 1.0 on conservation, and conservation alone would
    call it pristine.
    """
    if columns.empty:
        return None
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Scatter(x=columns["position"], y=columns["conservation"], mode="lines", name="Conservation", line=dict(color=_GREEN, width=1)),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(
            x=columns["position"],
            y=columns["gap_fraction"],
            mode="lines",
            name="Gap fraction",
            line=dict(color=_AMBER, width=1),
            fill="tozeroy",
            fillcolor="rgba(210,153,34,0.15)",
        ),
        secondary_y=True,
    )
    _theme(fig, height=320, title="Alignment profile")
    fig.update_xaxes(title_text="Alignment column")
    fig.update_yaxes(title_text="Conservation", range=[0, 1.05], secondary_y=False, gridcolor=_GRID)
    fig.update_yaxes(title_text="Gap fraction", range=[0, 1.05], secondary_y=True, showgrid=False)
    return fig


def branch_support_figure(support: "pd.DataFrame") -> Optional[go.Figure]:
    """Distribution of branch support, with the conventional 70% cut-off marked."""
    if support.empty:
        return None
    fig = go.Figure(go.Histogram(x=support["support"], marker_color=_BLUE, nbinsx=20))
    fig.add_vline(
        x=SUPPORT_THRESHOLD,
        line_dash="dash",
        line_color=_AMBER,
        annotation_text=f"{SUPPORT_THRESHOLD:.0f}% — clades below this are not considered supported",
        annotation_font_color=_AMBER,
    )
    _theme(fig, height=300, title="Branch support")
    fig.update_xaxes(title_text="Support (%)", range=[0, 105])
    fig.update_yaxes(title_text="Internal nodes")
    return fig


def unifrac_heatmap_figure(distances: "pd.DataFrame", weighted: bool = True) -> Optional[go.Figure]:
    """Between-sample UniFrac dissimilarity."""
    if distances.empty:
        return None
    samples = distances["sample"].tolist()
    values = distances[samples].values if all(name in distances for name in samples) else distances.iloc[:, 1:].values
    label = "Weighted UniFrac" if weighted else "Unweighted UniFrac"
    fig = go.Figure(
        go.Heatmap(z=values, x=samples, y=samples, colorscale="Viridis", colorbar=dict(title="distance"), hovertemplate="%{y} vs %{x}<br>%{z:.4f}<extra></extra>")
    )
    fig.update_layout(yaxis=dict(autorange="reversed"))
    subtitle = "abundance-weighted" if weighted else "presence/absence — rare lineages count as much as dominant ones"
    return _theme(fig, height=max(320, 16 * len(samples) + 140), title=f"{label} ({subtitle})")


def annotated_tree_figure(
    layout: dict,
    taxonomy: Optional[dict] = None,
    abundance: Optional[dict] = None,
) -> Optional[go.Figure]:
    """Phylogram with support-coloured nodes, taxonomy labels and abundance bars.

    The abundance bars share the tree's y-axis, so a tip and its bar are always on
    the same row.
    """
    if not layout or not layout.get("tips"):
        return None
    tips = layout["tips"]
    taxonomy = taxonomy or {}
    abundance = abundance or {}
    has_bars = bool(abundance)

    fig = (
        make_subplots(rows=1, cols=2, column_widths=[0.72, 0.28], shared_yaxes=True, horizontal_spacing=0.02)
        if has_bars
        else make_subplots(rows=1, cols=1)
    )

    line_x: list = []
    line_y: list = []
    for x_parent, x_child, y in layout.get("hlines", []):
        line_x += [x_parent, x_child, None]
        line_y += [y, y, None]
    for x, y0, y1 in layout.get("vlines", []):
        line_x += [x, x, None]
        line_y += [y0, y1, None]
    fig.add_trace(
        go.Scatter(x=line_x, y=line_y, mode="lines", line=dict(color=_TEAL, width=1), hoverinfo="skip", showlegend=False),
        row=1,
        col=1,
    )

    # Support values arrive on 0-1 (FastTree) or 0-100 (RAxML); decide from the tree.
    nodes = layout.get("nodes") or []
    if nodes:
        raw = [support for _, _, support in nodes]
        scale = 100.0 if max(raw) <= 1.0 else 1.0
        for supported, name, colour in ((True, f"support >= {SUPPORT_THRESHOLD:.0f}%", _GREEN), (False, f"support < {SUPPORT_THRESHOLD:.0f}%", _RED)):
            subset = [(x, y, support * scale) for x, y, support in nodes if (support * scale >= SUPPORT_THRESHOLD) == supported]
            if not subset:
                continue
            fig.add_trace(
                go.Scatter(
                    x=[point[0] for point in subset],
                    y=[point[1] for point in subset],
                    mode="markers",
                    marker=dict(size=7, color=colour, symbol="circle"),
                    name=name,
                    hovertext=[f"support {point[2]:.1f}%" for point in subset],
                    hovertemplate="%{hovertext}<extra></extra>",
                ),
                row=1,
                col=1,
            )

    labels = [f"{tip[2]} · {taxonomy[tip[2]]}" if taxonomy.get(tip[2]) else tip[2] for tip in tips]
    fig.add_trace(
        go.Scatter(
            x=[tip[0] for tip in tips],
            y=[tip[1] for tip in tips],
            mode="markers+text" if len(tips) <= 60 else "markers",
            text=labels if len(tips) <= 60 else None,
            textposition="middle right",
            textfont=dict(size=9),
            marker=dict(size=5, color=_BLUE),
            hovertext=labels,
            hovertemplate="%{hovertext}<extra></extra>",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    if has_bars:
        fig.add_trace(
            go.Bar(
                x=[abundance.get(tip[2], 0) for tip in tips],
                y=[tip[1] for tip in tips],
                orientation="h",
                marker_color=_PURPLE,
                name="Reads",
                hovertemplate="%{y}<br>%{x:,} reads<extra></extra>",
                showlegend=False,
            ),
            row=1,
            col=2,
        )
        fig.update_xaxes(title_text="Reads", row=1, col=2)

    _theme(fig, height=max(360, 16 * len(tips) + 140), title=f"Phylogenetic tree ({layout.get('n_tips', len(tips))} tips)")
    fig.update_xaxes(title_text="Branch length", row=1, col=1)
    fig.update_yaxes(showticklabels=False, showgrid=False)
    return fig




# ── Beta diversity / ordination ──────────────────────────────────────────

# Categorical palette for metadata groups. Distinguishable in both the dark theme and
# in greyscale print, which is where these figures usually end up.
_GROUP_COLORS = [_BLUE, _AMBER, _TEAL, _PURPLE, _GREEN, _RED, "#f778ba", "#8b98aa"]


def beta_ordination_figure(result: dict, groups: Optional[dict] = None) -> "go.Figure":
    """PC1 x PC2 scatter, labelled, with variance explained on the axes.

    The percentages are not decoration: without them the reader cannot tell whether an
    apparent separation lives on an axis carrying 60% of the variance or 3% of it, and
    the plot is uninterpretable.
    """
    coordinates = result["coordinates"]
    explained = result["explained"]
    samples = result["samples"]
    groups = groups or {}

    fig = go.Figure()
    labels = [groups.get(sample) or "no group" for sample in samples]
    for index, label in enumerate(dict.fromkeys(labels)):
        mask = [i for i, value in enumerate(labels) if value == label]
        fig.add_trace(
            go.Scatter(
                x=coordinates[mask, 0],
                y=coordinates[mask, 1],
                mode="markers+text",
                name=str(label),
                text=[samples[i] for i in mask],
                textposition="top center",
                textfont=dict(size=10, color="#c2ccd8"),
                marker=dict(size=13, color=_GROUP_COLORS[index % len(_GROUP_COLORS)], line=dict(width=1, color="#0b0f14")),
                hovertemplate="%{text}<br>PC1 %{x:.2f}<br>PC2 %{y:.2f}<extra></extra>",
            )
        )

    # Origin lines: in a CLR/PCoA space the centroid is the reference point, and reading
    # the quadrants is how you read the plot.
    fig.add_hline(y=0, line_width=1, line_dash="dot", line_color="rgba(255,255,255,0.15)")
    fig.add_vline(x=0, line_width=1, line_dash="dot", line_color="rgba(255,255,255,0.15)")

    name = DISTANCES.get(result["metric"], result["metric"])
    _theme(fig, height=430, title=f"{name} — {result['n_samples']} samples, {result['n_features']} features")
    fig.update_xaxes(title_text=f"PC1 ({explained[0] * 100:.1f}%)")
    fig.update_yaxes(title_text=f"PC2 ({explained[1] * 100:.1f}%)")
    # Equal aspect: distances in an ordination are only comparable if the axes are.
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    return fig


def scree_figure(result: dict) -> "go.Figure":
    """Variance explained per axis — the justification for how many axes to read."""
    explained = np.asarray(result["explained"], dtype=float)[:10] * 100
    axes = [f"PC{i}" for i in range(1, len(explained) + 1)]
    fig = go.Figure(
        go.Bar(
            x=axes,
            y=explained,
            marker_color=_BLUE,
            text=[f"{value:.1f}%" for value in explained],
            textposition="outside",
            hovertemplate="%{x}: %{y:.1f}%<extra></extra>",
        )
    )
    _theme(fig, height=260, title="Variance explained per axis")
    fig.update_yaxes(title_text="% of variance")
    return fig


def biplot_figure(result: dict, loadings: list[dict], groups: Optional[dict] = None) -> "go.Figure":
    """The ordination with feature vectors on top.

    This is what upgrades "the samples separate" into "the samples separate *because of
    this taxon*" — the vectors point toward the samples a feature is enriched in, and
    they are labelled by taxon rather than by OTU id wherever a name is known.
    """
    fig = beta_ordination_figure(result, groups)
    if not loadings:
        return fig

    coordinates = result["coordinates"]
    span = float(np.abs(coordinates).max()) or 1.0
    scale = span / (max(entry["magnitude"] for entry in loadings) or 1.0) * 0.8

    for entry in loadings:
        x, y = entry["pc1"] * scale, entry["pc2"] * scale
        fig.add_trace(
            go.Scatter(
                x=[0, x],
                y=[0, y],
                mode="lines+text",
                line=dict(color=_AMBER, width=1.5),
                text=["", entry["label"]],
                textposition="middle right",
                textfont=dict(size=9, color=_AMBER),
                hovertemplate=f"{entry['label']}<extra></extra>",
                showlegend=False,
            )
        )
    fig.update_layout(title=f"Biplot — top {len(loadings)} features")
    return fig


# ── Results summary ──────────────────────────────────────────────────────

# A step that loses more than this share of the reads is worth a red flag. The initial
# quality filter is exempt: discarding bad reads is what it is FOR.
FUNNEL_LOSS_ALERT = 0.20


def read_funnel_figure(steps: list[dict]) -> "go.Figure":
    """Every read from raw FASTQ to classified feature, in one picture.

    This is the plot that saves a user from visiting QC, Preprocess and Analysis in turn
    just to find out where their data went. Each bar carries the absolute count and the
    % of the original input, and any single step that drops more than 20% is coloured
    red — that is where to look first.
    """
    labels = [step["label"] for step in steps]
    values = [float(step["reads"]) for step in steps]
    baseline = values[0] if values else 1.0

    colors, annotations = [], []
    for index, value in enumerate(values):
        if index == 0:
            colors.append(_BLUE)
            annotations.append("")
            continue
        previous = values[index - 1] or 1.0
        loss = 1.0 - (value / previous)
        # The first filter is meant to lose reads; flagging it would cry wolf every run.
        alarming = loss > FUNNEL_LOSS_ALERT and index > 1
        colors.append(_RED if alarming else _TEAL)
        annotations.append(f"−{loss * 100:.0f}%" if loss > 0.005 else "")

    text = [
        f"{value:,.0f}<br>{value / baseline * 100:.0f}% of input" + (f"<br><b>{note}</b>" if note else "")
        for value, note in zip(values, annotations)
    ]

    fig = go.Figure(
        go.Bar(
            x=labels,
            y=values,
            marker_color=colors,
            text=text,
            textposition="outside",
            hovertemplate="%{x}<br>%{y:,.0f} reads<extra></extra>",
        )
    )
    _theme(fig, height=400, title="Read funnel — from primer removal to classified feature")
    fig.update_yaxes(title_text="Reads")
    # Headroom for the outside labels, which otherwise clip at the top.
    if values:
        fig.update_yaxes(range=[0, max(values) * 1.25])
    return fig


def stacked_composition_figure(frame: "pd.DataFrame", rank: str, top_n: int = 12) -> "go.Figure":
    """Stacked relative abundance per sample. Answers "what was found?" visually."""
    samples = [column for column in frame.columns if column != "taxon"]
    frame = frame.copy()
    frame["_total"] = frame[samples].sum(axis=1)
    frame = frame.sort_values("_total", ascending=False)

    # Beyond a dozen taxa the legend is unreadable and the colours stop being
    # distinguishable, so the tail is pooled honestly into "Outros".
    head = frame.head(top_n)
    tail = frame.iloc[top_n:]
    fig = go.Figure()
    for index, (_, row) in enumerate(head.iterrows()):
        fig.add_trace(
            go.Bar(
                x=samples,
                y=[row[sample] for sample in samples],
                name=str(row["taxon"])[:38],
                marker_color=_CATEGORICAL[index % (len(_CATEGORICAL) - 1)],
                hovertemplate="%{x}<br>" + str(row["taxon"]) + "<br>%{y:.1f}%<extra></extra>",
            )
        )
    if not tail.empty:
        fig.add_trace(
            go.Bar(
                x=samples,
                y=[tail[sample].sum() for sample in samples],
                name=f"Outros ({len(tail)} taxa)",
                marker_color="#8b98aa",
                hovertemplate="%{x}<br>Outros<br>%{y:.1f}%<extra></extra>",
            )
        )

    label = rank.replace("Consensus_", "")
    _theme(fig, height=420, title=f"Relative abundance per sample — {label}")
    fig.update_layout(barmode="stack", legend=dict(orientation="v", x=1.02, y=1, xanchor="left"))
    fig.update_yaxes(title_text="% of reads", range=[0, 100])
    return fig
