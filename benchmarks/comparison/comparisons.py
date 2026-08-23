"""How two outputs are compared. One copy, applied at evaluation time.

These used to live in four near-identical `_shared.py` modules and run on the GPU node,
which meant the definitions could drift apart and asking a new question of an old run cost
another GPU allocation. They are pure functions of stored arrays, so they belong here.

Each entry in `COMPARISONS` declares the arrays it needs at a comparison point. A criterion
named `<point>.<suffix>` is computed by looking `<suffix>` up here and loading those arrays
from the point. `reference` is always the CPU side, `candidate` the GPU side; several of
these are asymmetric, so the distinction matters.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from collections.abc import Callable

import numpy as np
from scipy import sparse
from sklearn.manifold import trustworthiness as _trustworthiness
from sklearn.metrics import accuracy_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors

RTOL = 1e-5
ATOL = 1e-8
KNN = 15


def dense(value) -> np.ndarray:
    if sparse.issparse(value):
        return value.toarray()
    return np.asarray(value)


def _floats(*values) -> tuple[np.ndarray, ...]:
    return tuple(dense(v).astype(float, copy=False) for v in values)


def _labels(value) -> np.ndarray:
    return np.asarray(dense(value)).astype(str).ravel()


# ---------------------------------------------------------------- numeric arrays


def allclose_excess(candidate, reference, *, rtol: float = RTOL, atol: float = ATOL) -> float:
    """Worst elementwise violation of `numpy.allclose`, as a fraction of its own envelope.

    `numpy.allclose` accepts when every element satisfies |a - b| <= atol + rtol * |b|, so
    dividing by that envelope gives one scale-free number: `<= 1` means the arrays are
    `allclose`. NaN propagates, because NaN is not close to anything and skipping it would
    make this weaker than the criterion it claims to be.
    """
    a, b = _floats(candidate, reference)
    return float(np.max(np.abs(a - b) / (atol + rtol * np.abs(b))))


def allclose_worst_magnitude(candidate, reference, *, rtol: float = RTOL, atol: float = ATOL) -> float:
    """Reference magnitude at the element that decides `allclose_excess`.

    This is what says which of the criterion's two terms bound it. Below `atol / rtol` the
    absolute floor dominates, and at float32 precision no implementation can satisfy it
    there; above, the disagreement is genuinely relative.
    """
    a, b = _floats(candidate, reference)
    ratio = (np.abs(a - b) / (atol + rtol * np.abs(b))).ravel()
    return float(np.abs(b.ravel()[int(np.nanargmax(ratio))]))


def allclose_violating_fraction(candidate, reference, *, rtol: float = RTOL, atol: float = ATOL) -> float:
    """Fraction of elements outside the envelope, not just the worst one.

    An elementwise criterion is decided by its tail, so the worst element alone cannot say
    whether two arrays disagree broadly or in a handful of places.
    """
    a, b = _floats(candidate, reference)
    return float(np.mean(np.abs(a - b) > atol + rtol * np.abs(b)))


def max_abs_error(candidate, reference) -> float:
    a, b = _floats(candidate, reference)
    return float(np.nanmax(np.abs(a - b)))


def mean_abs_error(candidate, reference) -> float:
    a, b = _floats(candidate, reference)
    return float(np.nanmean(np.abs(a - b)))


def max_rel_error(candidate, reference) -> float:
    """Largest absolute difference relative to the reference's own scale."""
    a, b = _floats(candidate, reference)
    scale = float(np.nanmax(np.abs(b)))
    return 0.0 if scale == 0.0 else float(np.nanmax(np.abs(a - b)) / scale)


def float32_ulp_at_max(candidate, reference) -> float:
    """One float32 ULP at the largest magnitude in the reference."""
    (b,) = _floats(reference)
    return float(np.spacing(np.float32(float(np.nanmax(np.abs(b))))))


def pearson_correlation(candidate, reference) -> float:
    """Correlation over the elements finite on both sides.

    Fewer than two such elements is not agreement, it is an absent measurement: NaN, which
    `evaluate.decide` fails. Returning a number there would score an all-NaN candidate as
    perfectly correlated with a finite reference. A constant pair is different — the
    correlation is undefined but the question "do they agree" still has an answer.
    """
    a, b = _floats(candidate, reference)
    a, b = a.ravel(), b.ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 2:
        return float("nan")
    if np.std(a[mask]) == 0 or np.std(b[mask]) == 0:
        return float(np.allclose(a[mask], b[mask]))
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def finite_mask_mean_abs_error(candidate, reference) -> float:
    """Mean absolute error over elements finite on both sides.

    p-values legitimately carry NaN where a test could not be run; those positions are
    checked separately by `nan_mask_agreement`.
    """
    a, b = _floats(candidate, reference)
    mask = np.isfinite(a) & np.isfinite(b)
    return float(np.mean(np.abs(a[mask] - b[mask]))) if mask.any() else 0.0


def nan_mask_agreement(candidate, reference) -> float:
    """Whether the two sides declare the same entries undefined."""
    a, b = _floats(candidate, reference)
    return float(np.mean(np.isnan(a) == np.isnan(b)))


def standard_deviation_max_abs_error(candidate, reference) -> float:
    a, b = _floats(candidate, reference)
    return float(np.nanmax(np.abs(a.std(axis=0) - b.std(axis=0))))


def relative_error_of_mean(candidate, reference) -> float:
    a, b = _floats(candidate, reference)
    return float(abs(a.mean() - b.mean()) / abs(b.mean()))


def relative_error_of_std(candidate, reference) -> float:
    a, b = _floats(candidate, reference)
    return float(abs(a.std() - b.std()) / abs(b.std()))


# ---------------------------------------------------------- components / embeddings


def _component_correlations(candidate, reference) -> np.ndarray:
    a, b = _floats(candidate, reference)
    n = min(a.shape[1], b.shape[1])
    return np.asarray([abs(pearson_correlation(a[:, i], b[:, i])) for i in range(n)])


def minimum_component_abs_correlation(candidate, reference) -> float:
    """Weakest per-component agreement. Sign is ignored: eigenvectors have no fixed sign."""
    return float(_component_correlations(candidate, reference).min())


def mean_component_abs_correlation(candidate, reference) -> float:
    return float(_component_correlations(candidate, reference).mean())


def _knn(values: np.ndarray, n_neighbors: int = KNN) -> np.ndarray:
    fitted = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(values)
    return fitted.kneighbors(return_distance=False)[:, 1:]


def cross_embedding_knn_overlap(candidate, reference, *, n_neighbors: int = KNN) -> float:
    """Mean per-cell overlap of the k nearest neighbours in two embeddings.

    Stochastic embeddings never match coordinate-wise, so agreement is measured as
    preservation of local neighbourhood structure.
    """
    a, b = _floats(candidate, reference)
    ka, kb = _knn(a, n_neighbors), _knn(b, n_neighbors)
    return float(np.mean([len(set(x) & set(y)) / n_neighbors for x, y in zip(ka, kb, strict=True)]))


def reference_trustworthiness(reference, reference_basis, *, n_neighbors: int = KNN) -> float:
    return float(_trustworthiness(dense(reference_basis), dense(reference), n_neighbors=n_neighbors))


def candidate_trustworthiness(candidate, candidate_basis, *, n_neighbors: int = KNN) -> float:
    return float(_trustworthiness(dense(candidate_basis), dense(candidate), n_neighbors=n_neighbors))


def trustworthiness_difference(candidate, candidate_basis, reference, reference_basis) -> float:
    return abs(
        candidate_trustworthiness(candidate, candidate_basis) - reference_trustworthiness(reference, reference_basis)
    )


def graph_jaccard(candidate, reference) -> float:
    """Mean per-row Jaccard overlap of two neighbour graphs, ignoring self-loops.

    Self is removed explicitly rather than by position, because the two implementations do
    not agree on whether the diagonal is stored.
    """
    a, b = candidate.tocsr(), reference.tocsr()
    scores = []
    for row in range(a.shape[0]):
        left = set(a.indices[a.indptr[row] : a.indptr[row + 1]]) - {row}
        right = set(b.indices[b.indptr[row] : b.indptr[row + 1]]) - {row}
        scores.append(len(left & right) / max(1, len(left | right)))
    return float(np.mean(scores))


# ------------------------------------------------------------------ labels / sets


def adjusted_rand_index(candidate, reference) -> float:
    return float(adjusted_rand_score(_labels(reference), _labels(candidate)))


def normalized_mutual_information(candidate, reference) -> float:
    return float(normalized_mutual_info_score(_labels(reference), _labels(candidate)))


def cluster_count_difference(candidate, reference) -> float:
    return float(abs(len(set(_labels(candidate))) - len(set(_labels(reference)))))


def exact_agreement(candidate, reference) -> float:
    """Fraction of positions where the two label vectors agree exactly."""
    a, b = _labels(candidate), _labels(reference)
    return float(np.mean(a == b))


def set_jaccard(candidate, reference) -> float:
    """Jaccard overlap of two name sets — gene selections, marker lists."""
    a, b = set(_labels(candidate)), set(_labels(reference))
    return len(a & b) / max(1, len(a | b))


def index_agreement(candidate, reference) -> float:
    """Whether two name vectors are identical in the same order."""
    a, b = _labels(candidate), _labels(reference)
    return float(a.shape == b.shape and bool(np.all(a == b)))


def label_nmi_reference(reference, reference_truth) -> float:
    return float(normalized_mutual_info_score(_labels(reference_truth), _labels(reference)))


def label_nmi_candidate(candidate, candidate_truth) -> float:
    return float(normalized_mutual_info_score(_labels(candidate_truth), _labels(candidate)))


def label_nmi_difference(reference, reference_truth, candidate, candidate_truth) -> float:
    return abs(label_nmi_candidate(candidate, candidate_truth) - label_nmi_reference(reference, reference_truth))


def reference_accuracy(reference, truth) -> float:
    return float(accuracy_score(_labels(truth), _labels(reference)))


def candidate_accuracy(candidate, truth) -> float:
    return float(accuracy_score(_labels(truth), _labels(candidate)))


def accuracy_difference(candidate, reference, truth) -> float:
    return abs(candidate_accuracy(candidate, truth) - reference_accuracy(reference, truth))


class Comparison(NamedTuple):
    inputs: tuple[str, ...]
    fn: Callable[..., float]


# Metric suffix -> the arrays it needs and how to combine them. A criterion named
# `<point>.<suffix>` is computed from the store whenever `<point>` carries those arrays;
# anything absent falls back to a scalar the measurement script recorded, so a committed
# archived run with no stored arrays still evaluates.
def _pair(fn):
    return Comparison(("candidate", "reference"), fn)


COMPARISONS: dict[str, Comparison] = {
    # numeric arrays
    "allclose_excess": _pair(allclose_excess),
    "pairwise_allclose_excess": _pair(allclose_excess),
    "score_allclose_excess": _pair(allclose_excess),
    "allclose_worst_magnitude": _pair(allclose_worst_magnitude),
    "allclose_violating_fraction": _pair(allclose_violating_fraction),
    "max_abs_error": _pair(max_abs_error),
    "mean_abs_error": _pair(finite_mask_mean_abs_error),
    "max_rel_error": _pair(max_rel_error),
    "float32_ulp_at_max": _pair(float32_ulp_at_max),
    "pearson_correlation": _pair(pearson_correlation),
    "pairwise_correlation": _pair(pearson_correlation),
    "score_correlation": _pair(pearson_correlation),
    "S_score_correlation": _pair(pearson_correlation),
    "G2M_score_correlation": _pair(pearson_correlation),
    "p_ko_correlation": _pair(pearson_correlation),
    "pca_correlation": _pair(pearson_correlation),
    "lda_abs_correlation": _pair(lambda c, r: abs(pearson_correlation(c, r))),
    "nan_mask_agreement": _pair(nan_mask_agreement),
    "standard_deviation_max_abs_error": _pair(standard_deviation_max_abs_error),
    "relative_error_of_mean": _pair(relative_error_of_mean),
    "relative_error_of_std": _pair(relative_error_of_std),
    # components and embeddings
    "minimum_component_abs_correlation": _pair(minimum_component_abs_correlation),
    "mean_component_abs_correlation": _pair(mean_component_abs_correlation),
    "cross_embedding_knn_overlap": _pair(cross_embedding_knn_overlap),
    "umap_knn_overlap": _pair(cross_embedding_knn_overlap),
    "graph_jaccard": _pair(graph_jaccard),
    "connectivity_jaccard": _pair(graph_jaccard),
    "cpu.trustworthiness": Comparison(("reference", "reference_basis"), reference_trustworthiness),
    "gpu.trustworthiness": Comparison(("candidate", "candidate_basis"), candidate_trustworthiness),
    "trustworthiness_difference": Comparison(
        ("candidate", "candidate_basis", "reference", "reference_basis"), trustworthiness_difference
    ),
    # labels and name sets
    "adjusted_rand_index": _pair(adjusted_rand_index),
    "normalized_mutual_information": _pair(normalized_mutual_information),
    "cluster_count_difference": _pair(cluster_count_difference),
    "exact_agreement": _pair(exact_agreement),
    "set_jaccard": _pair(set_jaccard),
    "index_agreement": _pair(index_agreement),
    "cpu_cell_type_nmi": Comparison(("reference", "reference_truth"), label_nmi_reference),
    "gpu_cell_type_nmi": Comparison(("candidate", "candidate_truth"), label_nmi_candidate),
    "cell_type_nmi_difference": Comparison(
        ("reference", "reference_truth", "candidate", "candidate_truth"), label_nmi_difference
    ),
    "cpu_accuracy": Comparison(("reference", "truth"), reference_accuracy),
    "gpu_accuracy": Comparison(("candidate", "truth"), candidate_accuracy),
    "accuracy_difference": Comparison(("candidate", "reference", "truth"), accuracy_difference),
}

# Extra evidence derived alongside a failing allclose criterion, so the reason a criterion
# failed never has to be anticipated at measurement time.
ALLCLOSE_EVIDENCE = (
    "allclose_worst_magnitude",
    "allclose_violating_fraction",
    "max_abs_error",
    "max_rel_error",
)

# Aggregations over metrics that already exist, for criteria stated over a family.
AGGREGATIONS: dict[str, tuple[str, Callable[[list[float]], float]]] = {
    "mean_set_jaccard": ("set_jaccard", lambda values: float(np.mean(values))),
    "minimum_set_jaccard": ("set_jaccard", lambda values: float(np.min(values))),
}
