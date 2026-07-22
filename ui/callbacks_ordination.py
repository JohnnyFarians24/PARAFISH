"""Ordination callbacks, registered once per mount point.

One implementation drives Analysis (full controls), Results (compact) and Phylo
(UniFrac pre-selected). ``register(app, prefix=...)`` wires a mount.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from dash import Input, Output, State, dcc, html, no_update
from dash.exceptions import PreventUpdate

from analysis.ordination import (
    LOW_CONFIDENCE_SAMPLES,
    NEEDS_TREE,
    OrdinationError,
    ordinate,
    permanova,
    top_loadings,
    write_coordinates,
    write_distance_matrix,
)
from ui.components.ordination_panel import guard_message, low_confidence_warning
from ui.components.table import data_table

_GRAPH_CONFIG = {"displaylogo": False, "responsive": True}

# Taxa that are not the biology under study. The host fish dominates an 18S library by
# orders of magnitude, so leaving it in guarantees PC1 is "how much host DNA did this
# sample happen to have" — which is a library-prep artefact, not an ecological gradient.
_HOST_MARKERS = ("chordata", "craniata", "actinopteri", "vertebrata", "mammalia")


def _run_dir(run_config_data, project_data) -> Optional[Path]:
    from core.run_registry import run_dir_for

    return run_dir_for(run_config_data, project_data)


def _matrix_path(run_dir: Path, feature_set: str = "otu") -> Optional[Path]:
    if not run_dir:
        return None
    for candidate in (
        run_dir / "analysis" / "otu_matrix.csv",
        run_dir / "analysis" / "pooled" / "otu_table.tsv",
        run_dir / "analysis" / "dada2" / "asv_abundance.csv",
    ):
        if candidate.exists():
            if feature_set == "asv" and "dada2" not in candidate.parts:
                continue
            if feature_set == "otu" and "dada2" in candidate.parts:
                continue
            return candidate
    return None


def _taxonomy_map(run_dir: Path, feature_set: str = "otu") -> dict[str, str]:
    """feature id -> best available taxon name, for labelling biplot vectors."""
    import csv

    if not run_dir:
        return {}
    path = run_dir / "databases" / feature_set / "taxonomy_consensus.csv"
    if not path.exists():
        path = run_dir / "databases" / "taxonomy_consensus.csv"
    if not path.exists():
        return {}

    names: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            feature = (row.get("ASV_ID") or "").strip()
            if not feature:
                continue
            for column in ("Consensus_Species", "Consensus_Genus", "Consensus_Family", "Consensus_Phylum"):
                value = (row.get(column) or "").strip()
                if value:
                    names[feature] = value
                    break
    return names


def _host_features(run_dir: Path, feature_set: str = "otu") -> set[str]:
    """Features assigned to the host lineage, or to nothing at all."""
    import csv

    if not run_dir:
        return set()
    path = run_dir / "databases" / feature_set / "taxonomy_consensus.csv"
    if not path.exists():
        return set()

    excluded: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            feature = (row.get("ASV_ID") or "").strip()
            if not feature:
                continue
            lineage = " ".join(
                (row.get(column) or "") for column in ("Consensus_Phylum", "Consensus_Class", "Consensus_Order")
            ).lower()
            classified = any((row.get(column) or "").strip() for column in ("Consensus_Genus", "Consensus_Family"))
            if any(marker in lineage for marker in _HOST_MARKERS) or not classified:
                excluded.add(feature)
    return excluded


def _sample_groups(samples_data) -> dict[str, str]:
    groups = {}
    for record in (samples_data or {}).get("samples", []) or []:
        sample_id = record.get("internal_id")
        if not sample_id:
            continue
        if record.get("is_control"):
            groups[sample_id] = "controlo"
        elif record.get("group"):
            groups[sample_id] = str(record["group"])
    return groups


def compute(run_config_data, project_data, samples_data, metric, prevalence, exclude_host, feature_set="otu"):
    """Run the ordination for the active run. Returns ``(result, groups, taxonomy, error)``."""
    from analysis.analysis_metrics import load_matrix

    run_dir = _run_dir(run_config_data, project_data)
    matrix_path = _matrix_path(run_dir, feature_set) if run_dir else None
    if not matrix_path:
        return None, {}, {}, "Sem tabela de abundância. Corre a etapa Analysis primeiro."

    features, samples, counts = load_matrix(matrix_path)
    if not features or not samples:
        return None, {}, {}, "A tabela de abundância está vazia."

    if exclude_host:
        drop = _host_features(run_dir, feature_set)
        keep = [index for index, feature in enumerate(features) if feature not in drop]
        if len(keep) >= 3:
            counts = counts[:, keep]
            features = [features[index] for index in keep]

    tree_path = None
    if metric in NEEDS_TREE:
        from ui.feature_sets import phylo_tree_path

        resolved = phylo_tree_path(run_dir, feature_set)
        if resolved:
            tree_path = Path(resolved)

    try:
        result = ordinate(
            counts,
            samples,
            features,
            metric=metric or "aitchison",
            min_prevalence=int(prevalence or 2),
            tree_path=tree_path,
        )
    except OrdinationError as exc:
        return None, {}, {}, str(exc)

    return result, _sample_groups(samples_data), _taxonomy_map(run_dir, feature_set), None


def register(app: Any, prefix: str = "ord-analysis", *, compact: bool = False) -> None:
    """Wire one ordination mount."""

    outputs = [Output(f"{prefix}-plot", "children")]
    if not compact:
        outputs += [
            Output(f"{prefix}-biplot", "children"),
            Output(f"{prefix}-scree", "children"),
            Output(f"{prefix}-permanova", "children"),
        ]

    inputs = [
        Input("store-analysis-outputs", "data"),
        Input(f"{prefix}-metric", "value") if not compact else Input("store-run-config", "data"),
    ]
    if not compact:
        inputs += [
            Input(f"{prefix}-prevalence", "value"),
            Input(f"{prefix}-exclude-host", "value"),
        ]

    states = [State("store-run-config", "data"), State("store-project", "data"), State("store-samples", "data")]

    @app.callback(*outputs, *inputs, *states, prevent_initial_call=False)
    def _render(*args):
        if compact:
            _analysis, _cfg = args[0], args[1]
            metric, prevalence, exclude = "aitchison", 2, ["exclude"]
            run_config_data, project_data, samples_data = args[2:5]
        else:
            _analysis, metric, prevalence, exclude = args[0:4]
            run_config_data, project_data, samples_data = args[4:7]

        result, groups, taxonomy, error = compute(
            run_config_data, project_data, samples_data, metric, prevalence, bool(exclude)
        )

        if error:
            # Never a plot when a plot would mislead. The message says why, in full.
            blank = guard_message(0, error)
            return (blank,) if compact else (blank, blank, blank, blank)

        from ui.plots import beta_ordination_figure, biplot_figure, scree_figure

        children = []
        if result["low_confidence"]:
            children.append(low_confidence_warning(result["n_samples"]))
        children.append(dcc.Graph(figure=beta_ordination_figure(result, groups), config=_GRAPH_CONFIG))
        plot = html.Div(children)
        if compact:
            return (plot,)

        loadings = top_loadings(result, taxonomy, top_n=10)
        if loadings:
            biplot = html.Div(
                [
                    dcc.Graph(figure=biplot_figure(result, loadings, groups), config=_GRAPH_CONFIG),
                    data_table(
                        data=[
                            {
                                "taxon": entry["label"],
                                "feature": entry["feature"],
                                "PC1": round(entry["pc1"], 3),
                                "PC2": round(entry["pc2"], 3),
                            }
                            for entry in loadings
                        ],
                        columns=[{"name": n, "id": i} for n, i in [("Taxon", "taxon"), ("Feature", "feature"), ("PC1", "PC1"), ("PC2", "PC2")]],
                    ),
                ]
            )
        else:
            biplot = html.Div("O biplot só está disponível para a PCA de Aitchison (as PCoA não têm cargas por feature).", className="empty-state-sub")

        scree = dcc.Graph(figure=scree_figure(result), config=_GRAPH_CONFIG)

        stats = permanova(result, groups) if groups else None
        if stats is None:
            perm = html.Div(
                "A PERMANOVA precisa de pelo menos 2 grupos com 3 ou mais réplicas cada. "
                "Define grupos na tabela de amostras para a activar.",
                className="empty-state-sub",
            )
        else:
            floor = stats.get("p_floor")
            note = ""
            if floor and stats["p_value"] <= floor * 1.5:
                # A permutation test cannot return a p below 1/(number of relabellings).
                note = (
                    f" Nota: com estes tamanhos de grupo o p mínimo possível é {floor:.2f}, "
                    f"por isso este valor está no chão do teste e não deve ser lido como ausência de efeito."
                )
            perm = html.Div(
                [
                    html.Div(
                        [
                            html.Div(f"pseudo-F = {stats['pseudo_f']:.2f}", className="status-pill"),
                            html.Div(f"R² = {stats['r_squared']:.3f}", className="status-pill"),
                            html.Div(f"p = {stats['p_value']:.3f}", className="status-pill"),
                        ],
                        className="chip-list",
                    ),
                    html.Div(
                        f"{stats['permutations']} permutações; grupos: "
                        + ", ".join(f"{name} (n={n})" for name, n in stats["groups"].items())
                        + note,
                        className="mini-note",
                    ),
                ]
            )

        return plot, biplot, scree, perm

    if compact:
        return

    @app.callback(
        Output(f"{prefix}-download-coords-data", "data"),
        Input(f"{prefix}-download-coords", "n_clicks"),
        State(f"{prefix}-metric", "value"),
        State(f"{prefix}-prevalence", "value"),
        State(f"{prefix}-exclude-host", "value"),
        State("store-run-config", "data"),
        State("store-project", "data"),
        State("store-samples", "data"),
        prevent_initial_call=True,
    )
    def _download_coords(n_clicks, metric, prevalence, exclude, run_config_data, project_data, samples_data):
        if not n_clicks:
            raise PreventUpdate
        result, _groups, _tax, error = compute(run_config_data, project_data, samples_data, metric, prevalence, bool(exclude))
        if error or result is None:
            raise PreventUpdate
        import tempfile

        path = Path(tempfile.mkdtemp()) / "ordination_coordinates.tsv"
        write_coordinates(result, path)
        return dcc.send_file(str(path))

    @app.callback(
        Output(f"{prefix}-download-dist-data", "data"),
        Input(f"{prefix}-download-dist", "n_clicks"),
        State(f"{prefix}-metric", "value"),
        State(f"{prefix}-prevalence", "value"),
        State(f"{prefix}-exclude-host", "value"),
        State("store-run-config", "data"),
        State("store-project", "data"),
        State("store-samples", "data"),
        prevent_initial_call=True,
    )
    def _download_dist(n_clicks, metric, prevalence, exclude, run_config_data, project_data, samples_data):
        if not n_clicks:
            raise PreventUpdate
        result, _groups, _tax, error = compute(run_config_data, project_data, samples_data, metric, prevalence, bool(exclude))
        if error or result is None:
            raise PreventUpdate
        import tempfile

        path = Path(tempfile.mkdtemp()) / "ordination_distances.tsv"
        write_distance_matrix(result, path)
        return dcc.send_file(str(path))

    @app.callback(
        Output(f"{prefix}-colour", "options"),
        Output(f"{prefix}-colour", "value"),
        Input("store-samples", "data"),
    )
    def _colour_options(samples_data):
        groups = _sample_groups(samples_data)
        if groups:
            return [{"label": "Grupo de metadados", "value": "group"}, {"label": "Amostra", "value": "sample"}], "group"
        # No metadata to colour by: say so instead of offering a dead dropdown.
        return [{"label": "Amostra (sem grupos definidos)", "value": "sample"}], "sample"
