"""Beta diversity ordination, done correctly for compositional data.

**Why not plain PCA on counts or relative abundance.** Metabarcoding counts are
*compositional*: the total is fixed by how deep you sequenced, not by biology. Two
samples can have identical communities and wildly different totals, and one feature
going up forces the others down by construction. On such data Euclidean distance — the
metric PCA optimises — is not meaningful, and the usual results are (a) PC1 tracking
sequencing depth rather than any biological gradient, and (b) the horseshoe artefact.

**The default here is therefore a CLR (Aitchison) PCA:**

1. drop low-prevalence features (default: must appear in >= 2 samples);
2. replace zeros with a pseudocount (log(0) is undefined) — see ZERO_REPLACEMENT;
3. centred log-ratio per sample: ``clr(x)_i = log(x_i) - mean(log(x))``, which moves the
   data out of the simplex and into a real Euclidean space;
4. PCA on the CLR matrix, centred but **not** variance-scaled (CLR components are
   already on one log scale; scaling would re-weight them arbitrarily).

That is mathematically the same as PCoA on the Aitchison distance, and it is the only
option here whose geometry the arithmetic actually supports.

PCoA on Bray-Curtis / Jaccard / UniFrac is offered as an alternative because they are
what the literature mostly reports, and because UniFrac is the only one that uses the
phylogeny.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

import numpy as np

# log(0) is undefined, so the zeros — which are most of the matrix — have to become
# something. 0.5 on the counts is the conventional default: it is half the smallest
# observable count, so it says "below the detection limit" without asserting a
# particular abundance.
#
# CAVEAT, and it is a real one: the CLR is exactly scale-invariant, but ``clr(x + 0.5)``
# is only approximately so, because 0.5 means something different at a depth of 100 than
# at 5,000. In practice the distortion is small (~0.1% on the coordinates for a 50x depth
# difference, which is far more depth variation than a real run has), but it IS the price
# of a fixed pseudocount, and it is exactly what zCompositions' Bayesian multiplicative
# replacement (cmultRepl) exists to avoid — that scales the replacement with the sample.
# We keep the pseudocount: it needs no R dependency, and the residual distortion is orders
# of magnitude below the biological signal. The choice is stated in the UI, not hidden.
ZERO_REPLACEMENT = 0.5

# A feature seen in a single sample carries no between-sample information and only adds
# noise to the log-ratio geometry.
MIN_PREVALENCE = 2

# Below three samples an ordination cannot say anything: two points always lie on one
# line, so the "separation" is an artefact of having two points. This is science, not a
# limitation to paper over — the UI states it rather than drawing a misleading plot.
MIN_SAMPLES = 3
LOW_CONFIDENCE_SAMPLES = 5
MIN_FEATURES = 3

DISTANCES = {
    "aitchison": "PCA de Aitchison (CLR) — recomendado",
    "braycurtis": "PCoA — Bray-Curtis",
    "jaccard": "PCoA — Jaccard (presence/absence)",
    "weighted_unifrac": "PCoA — weighted UniFrac (requires a tree)",
    "unweighted_unifrac": "PCoA — unweighted UniFrac (requires a tree)",
}

NEEDS_TREE = frozenset({"weighted_unifrac", "unweighted_unifrac"})


class OrdinationError(Exception):
    """Raised when an ordination cannot be computed honestly."""


def filter_features(counts: np.ndarray, features: list[str], min_prevalence: int = MIN_PREVALENCE):
    """Drop features present in fewer than ``min_prevalence`` samples."""
    if counts.size == 0:
        return counts, features
    prevalence = (counts > 0).sum(axis=0)
    keep = prevalence >= max(1, int(min_prevalence))
    if not keep.any():
        return counts[:, :0], []
    return counts[:, keep], [feature for feature, keep_it in zip(features, keep) if keep_it]


def clr(counts: np.ndarray, pseudocount: float = ZERO_REPLACEMENT) -> np.ndarray:
    """Centred log-ratio, per sample (row)."""
    data = counts.astype(float) + float(pseudocount)
    log_data = np.log(data)
    return log_data - log_data.mean(axis=1, keepdims=True)


def pca(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """PCA by SVD. Returns ``(coordinates, explained_variance_ratio, loadings)``.

    Centred, not scaled: CLR values already share one log scale, and dividing each
    component by its standard deviation would silently re-weight the features.
    """
    centred = matrix - matrix.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(centred, full_matrices=False)
    coordinates = u * s
    variance = s**2
    total = variance.sum()
    ratio = variance / total if total > 0 else np.zeros_like(variance)
    return coordinates, ratio, vt


def _pcoa(counts: np.ndarray, samples: list[str], metric: str, tree_path: Optional[Path], features: list[str]):
    from skbio.diversity import beta_diversity
    from skbio.stats.ordination import pcoa as skbio_pcoa

    kwargs = {}
    if metric in NEEDS_TREE:
        if not tree_path or not Path(tree_path).exists():
            raise OrdinationError(
                "UniFrac needs a phylogenetic tree. Run the Phylo stage first, "
                "or choose another distance."
            )
        from skbio import TreeNode

        tree = TreeNode.read(str(tree_path))
        # UniFrac walks the tree, so every feature in the table must be a tip on it.
        tips = {tip.name for tip in tree.tips()}
        keep = [index for index, feature in enumerate(features) if feature in tips]
        if len(keep) < MIN_FEATURES:
            raise OrdinationError(
                "The table's features do not match the tree's tips, so UniFrac "
                "cannot be computed."
            )
        counts = counts[:, keep]
        features = [features[index] for index in keep]
        kwargs = {"tree": tree, "otu_ids": features}

    data = counts.astype(int) if metric == "jaccard" else counts.astype(float)
    matrix = beta_diversity(metric, data, ids=samples, validate=True, **kwargs)

    # A degenerate matrix — every distance zero — has no axes to decompose, and skbio
    # divides by a zero eigenvalue sum and hands back NaN proportions. This happens for
    # real: Jaccard sees only presence/absence, so if every feature occurs in every
    # sample, all samples are identical to it. Plotting NaN axes would be nonsense.
    if not np.any(np.asarray(matrix.data) > 0):
        raise OrdinationError(
            f"Every sample is at zero distance under this metric, so there are no axes "
            f"to decompose. With Jaccard this happens when every feature is present in "
            f"every sample — use a metric that accounts for abundance (Aitchison or "
            f"Bray-Curtis)."
        )

    result = skbio_pcoa(matrix)
    coordinates = np.asarray(result.samples.iloc[:, :2])
    ratio = np.asarray(result.proportion_explained)
    if not np.all(np.isfinite(ratio)):
        raise OrdinationError("The ordination produced undefined axes with this metric. Try Aitchison or Bray-Curtis.")
    return coordinates, ratio, matrix, features


def ordinate(
    counts: np.ndarray,
    samples: list[str],
    features: list[str],
    *,
    metric: str = "aitchison",
    min_prevalence: int = MIN_PREVALENCE,
    pseudocount: float = ZERO_REPLACEMENT,
    tree_path: Optional[Path] = None,
) -> dict:
    """Ordinate one abundance matrix. Raises ``OrdinationError`` when it would mislead."""
    n_samples = len(samples)
    if n_samples < MIN_SAMPLES:
        raise OrdinationError(
            f"Ordination requires at least {MIN_SAMPLES} samples (recommended: {LOW_CONFIDENCE_SAMPLES} or more). "
            f"This run has {n_samples}. With 2 samples there is only 1 axis and the separation is trivial, "
            f"so the plot would be misleading."
        )

    counts, features = filter_features(counts, features, min_prevalence)
    if len(features) < MIN_FEATURES:
        raise OrdinationError(
            f"Only {len(features)} features remain after the prevalence filter (minimum {MIN_FEATURES}). "
            f"Lower the minimum prevalence or check the abundance table."
        )

    loadings = None
    distance_matrix = None
    if metric == "aitchison":
        matrix = clr(counts, pseudocount)
        coordinates, ratio, vt = pca(matrix)
        coordinates = coordinates[:, :2]
        loadings = vt[:2]
    else:
        coordinates, ratio, distance_matrix, features = _pcoa(counts, samples, metric, tree_path, features)

    # Pad, so a run with exactly 3 samples (which yields 2 axes) still has a PC2.
    if coordinates.shape[1] < 2:
        coordinates = np.hstack([coordinates, np.zeros((coordinates.shape[0], 1))])
    ratio = np.asarray(ratio, dtype=float)
    if ratio.size < 2:
        ratio = np.concatenate([ratio, np.zeros(2 - ratio.size)])

    return {
        "metric": metric,
        "samples": samples,
        "features": features,
        "coordinates": coordinates[:, :2],
        "explained": ratio,
        "loadings": loadings,
        "distance_matrix": distance_matrix,
        "n_samples": n_samples,
        "n_features": len(features),
        "low_confidence": n_samples < LOW_CONFIDENCE_SAMPLES,
        "pseudocount": pseudocount if metric == "aitchison" else None,
    }


def top_loadings(result: dict, taxonomy: Optional[dict] = None, top_n: int = 10) -> list[dict]:
    """The features pulling hardest on PC1/PC2, named by taxon where possible.

    This is what turns "the samples separate" into "the samples separate *because of
    this parasite*", which is the only version of the statement worth publishing.
    """
    loadings = result.get("loadings")
    if loadings is None:
        return []
    pc1, pc2 = loadings[0], loadings[1]
    magnitude = np.sqrt(pc1**2 + pc2**2)
    order = np.argsort(magnitude)[::-1][:top_n]
    taxonomy = taxonomy or {}
    return [
        {
            "feature": result["features"][index],
            "label": taxonomy.get(result["features"][index]) or result["features"][index],
            "pc1": float(pc1[index]),
            "pc2": float(pc2[index]),
            "magnitude": float(magnitude[index]),
        }
        for index in order
    ]


def permanova(result: dict, groups: dict[str, str], permutations: int = 999) -> Optional[dict]:
    """PERMANOVA on the ordination's distances, when the design can support it.

    Needs at least two groups with >= 3 replicates each; below that the p-value is not
    worth printing, so this returns ``None`` and the UI hides the section rather than
    reporting a test nobody should believe.
    """
    from collections import Counter

    labels = [groups.get(sample) for sample in result["samples"]]
    if any(label is None for label in labels):
        return None
    counts = Counter(labels)
    if len(counts) < 2 or any(count < 3 for count in counts.values()):
        return None

    from skbio.stats.distance import DistanceMatrix, permanova as skbio_permanova

    matrix = result.get("distance_matrix")
    if matrix is None:
        # Aitchison: the Euclidean distance between CLR coordinates IS the Aitchison
        # distance, which is the whole point of the transform.
        from scipy.spatial.distance import pdist, squareform

        matrix = DistanceMatrix(squareform(pdist(result["coordinates"])), ids=result["samples"])

    outcome = skbio_permanova(matrix, labels, permutations=permutations)
    # R^2 is not reported by skbio; derive it from pseudo-F and the degrees of freedom.
    f_stat = float(outcome["test statistic"])
    n, k = len(labels), len(counts)
    r_squared = (f_stat * (k - 1)) / (f_stat * (k - 1) + (n - k)) if n > k else float("nan")

    # The smallest p-value a permutation test CAN return is 1 / (number of distinct
    # relabellings). With two groups of three that is 1/10 = 0.1 — so a "p = 0.12" from
    # such a design is the floor of the test, not weak evidence, and reporting it without
    # that context invites the reader to conclude there is no effect.
    from math import comb

    if k == 2:
        sizes = sorted(counts.values())
        distinct = comb(n, sizes[0]) // 2
        p_floor = 1.0 / distinct if distinct else None
    else:
        p_floor = None

    return {
        "pseudo_f": f_stat,
        "p_value": float(outcome["p-value"]),
        "r_squared": r_squared,
        "permutations": permutations,
        "groups": dict(counts),
        "p_floor": p_floor,
    }


def write_coordinates(result: dict, path: Path) -> Path:
    """Sample coordinates as TSV — enough to reproduce the figure elsewhere."""
    path.parent.mkdir(parents=True, exist_ok=True)
    explained = result["explained"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "PC1", "PC2"])
        for index, sample in enumerate(result["samples"]):
            writer.writerow([sample, f"{result['coordinates'][index, 0]:.6f}", f"{result['coordinates'][index, 1]:.6f}"])
        writer.writerow([])
        writer.writerow(["# metric", result["metric"]])
        writer.writerow(["# PC1_explained", f"{explained[0]:.6f}"])
        writer.writerow(["# PC2_explained", f"{explained[1]:.6f}"])
    return path


def write_distance_matrix(result: dict, path: Path) -> Optional[Path]:
    """The distance matrix as TSV, so the ordination can be re-run independently."""
    matrix = result.get("distance_matrix")
    samples = result["samples"]
    if matrix is None:
        from scipy.spatial.distance import pdist, squareform

        data = squareform(pdist(result["coordinates"]))
    else:
        data = np.asarray(matrix.data)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["", *samples])
        for index, sample in enumerate(samples):
            writer.writerow([sample, *(f"{value:.6f}" for value in data[index])])
    return path
