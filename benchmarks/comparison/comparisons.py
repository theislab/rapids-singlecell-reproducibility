from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from collections.abc import Callable

import numpy as np
from scipy import sparse
from sklearn.manifold import trustworthiness as _trustworthiness
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors

RTOL = 1e-5
ATOL = 1e-8  # numpy.allclose default
FLOAT32_ATOL = 1e-6  # the default ATOL is below float32 resolution near zero
KNN = 15


def dense(value) -> np.ndarray:
    if sparse.issparse(value):
        return value.toarray()
    return np.asarray(value)


def _floats(*values) -> tuple[np.ndarray, ...]:
    arrays = tuple(dense(v).astype(float, copy=False) for v in values)
    if len({a.shape for a in arrays}) > 1:
        raise ValueError(f"shapes differ: {[a.shape for a in arrays]}")
    return arrays


def _labels(value) -> np.ndarray:
    return np.asarray(dense(value)).astype(str).ravel()


# ---------------------------------------------------------------- numeric arrays


def allclose_excess(candidate, reference, *, atol: float = ATOL) -> float:
    a, b = _floats(candidate, reference)
    return float(np.max(np.abs(a - b) / (atol + RTOL * np.abs(b))))


def pearson_correlation(candidate, reference) -> float:
    a, b = _floats(candidate, reference)
    a, b = a.ravel(), b.ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 2:
        return float("nan")
    if np.std(a[mask]) == 0 or np.std(b[mask]) == 0:
        return float(np.allclose(a[mask], b[mask]))
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def finite_mask_mean_abs_error(candidate, reference) -> float:
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


def relative_l2_max(candidate, reference) -> float:
    """Worst per-component relative L2 error, ||a - b|| / ||b||."""
    a, b = _floats(candidate, reference)
    denominator = np.linalg.norm(b, axis=0)
    numerator = np.linalg.norm(a - b, axis=0)
    # A zero-norm reference component is degenerate: any difference is infinitely relative,
    # and calling that agreement would hide it.
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(denominator > 0, numerator / denominator, np.where(numerator > 0, np.inf, 0.0))
    return float(np.nanmax(ratio))


def relative_error_of_mean(candidate, reference) -> float:
    a, b = _floats(candidate, reference)
    return float(abs(a.mean() - b.mean()) / abs(b.mean()))


def relative_error_of_std(candidate, reference) -> float:
    a, b = _floats(candidate, reference)
    return float(abs(a.std() - b.std()) / abs(b.std()))


# ---------------------------------------------------------- components / embeddings


def _component_correlations(candidate, reference) -> np.ndarray:
    a, b = _floats(candidate, reference)
    return np.asarray([abs(pearson_correlation(a[:, i], b[:, i])) for i in range(b.shape[1])])


def abs_allclose_excess(candidate, reference, *, atol: float = ATOL) -> float:
    """allclose_excess on magnitudes, because eigenvector signs are arbitrary."""
    a, b = _floats(candidate, reference)
    return allclose_excess(np.abs(a), np.abs(b), atol=atol)


def minimum_component_abs_correlation(candidate, reference) -> float:
    """Weakest per-component agreement. Sign is ignored: eigenvectors have no fixed sign."""
    return float(_component_correlations(candidate, reference).min())


def mean_component_abs_correlation(candidate, reference) -> float:
    return float(_component_correlations(candidate, reference).mean())


def _knn(values: np.ndarray) -> np.ndarray:
    fitted = NearestNeighbors(n_neighbors=KNN + 1).fit(values)
    return fitted.kneighbors(return_distance=False)[:, 1:]


def cross_embedding_knn_overlap(candidate, reference) -> float:
    """Mean per-cell overlap of the k nearest neighbours in two embeddings."""
    a, b = _floats(candidate, reference)
    return float(np.mean([len(set(x) & set(y)) / KNN for x, y in zip(_knn(a), _knn(b), strict=True)]))


def trustworthiness(embedding, basis) -> float:
    return float(_trustworthiness(dense(basis), dense(embedding), n_neighbors=KNN))


def _row_sets(candidate, reference):
    """Neighbour sets per row of two graphs, self-loops removed."""
    a, b = candidate.tocsr(), reference.tocsr()
    for row in range(a.shape[0]):
        yield (
            set(a.indices[a.indptr[row] : a.indptr[row + 1]]) - {row},
            set(b.indices[b.indptr[row] : b.indptr[row + 1]]) - {row},
        )


def graph_jaccard(candidate, reference) -> float:
    return float(np.mean([len(x & y) / len(x | y) if x | y else 1.0 for x, y in _row_sets(candidate, reference)]))


# ------------------------------------------------------------------ labels / sets


def adjusted_rand_index(candidate, reference) -> float:
    return float(adjusted_rand_score(_labels(reference), _labels(candidate)))


def normalized_mutual_information(candidate, reference) -> float:
    return float(normalized_mutual_info_score(_labels(reference), _labels(candidate)))


def exact_agreement(candidate, reference) -> float:
    """Fraction of positions where the two vectors agree exactly; numbers compare as numbers."""
    a, b = dense(candidate).ravel(), dense(reference).ravel()
    if a.dtype.kind not in "biuf" or b.dtype.kind not in "biuf":
        a, b = _labels(a), _labels(b)
    return float(np.mean(a == b))


def set_jaccard(candidate, reference) -> float:
    """Jaccard overlap of two name sets — gene selections, marker lists."""
    a, b = set(_labels(candidate)), set(_labels(reference))
    return len(a & b) / len(a | b) if a | b else 1.0


def graph_exact_agreement(candidate, reference) -> float:
    return float(np.mean([x == y for x, y in _row_sets(candidate, reference)]))


def index_agreement(candidate, reference) -> float:
    """Whether two name vectors are identical in the same order."""
    a, b = _labels(candidate), _labels(reference)
    return float(a.shape == b.shape and bool(np.all(a == b)))


def label_nmi_difference(candidate, reference, truth) -> float:
    return abs(normalized_mutual_information(candidate, truth) - normalized_mutual_information(reference, truth))


def accuracy_difference(candidate, reference, truth) -> float:
    return abs(exact_agreement(candidate, truth) - exact_agreement(reference, truth))


class Comparison(NamedTuple):
    inputs: tuple[str, ...]
    fn: Callable[..., float]


# Metric suffix -> the arrays it needs at a comparison point, and the function on them.
def _pair(fn):
    return Comparison(("candidate", "reference"), fn)


COMPARISONS: dict[str, Comparison] = {
    # numeric arrays
    "allclose_excess": _pair(allclose_excess),
    "mean_abs_error": _pair(finite_mask_mean_abs_error),
    "pearson_correlation": _pair(pearson_correlation),
    "lda_abs_correlation": _pair(lambda c, r: abs(pearson_correlation(c, r))),
    "nan_mask_agreement": _pair(nan_mask_agreement),
    "standard_deviation_max_abs_error": _pair(standard_deviation_max_abs_error),
    "relative_l2_max": _pair(relative_l2_max),
    "abs_allclose_excess": _pair(abs_allclose_excess),
    "float32_allclose_excess": _pair(lambda c, r: allclose_excess(c, r, atol=FLOAT32_ATOL)),
    "abs_float32_allclose_excess": _pair(lambda c, r: abs_allclose_excess(c, r, atol=FLOAT32_ATOL)),
    "relative_error_of_mean": _pair(relative_error_of_mean),
    "relative_error_of_std": _pair(relative_error_of_std),
    # components and embeddings
    "minimum_component_abs_correlation": _pair(minimum_component_abs_correlation),
    "mean_component_abs_correlation": _pair(mean_component_abs_correlation),
    "cross_embedding_knn_overlap": _pair(cross_embedding_knn_overlap),
    "graph_jaccard": _pair(graph_jaccard),
    "graph_exact_agreement": _pair(graph_exact_agreement),
    "cpu.trustworthiness": Comparison(("reference", "basis"), trustworthiness),
    "gpu.trustworthiness": Comparison(("candidate", "basis"), trustworthiness),
    # labels and name sets
    "adjusted_rand_index": _pair(adjusted_rand_index),
    "normalized_mutual_information": _pair(normalized_mutual_information),
    "exact_agreement": _pair(exact_agreement),
    "set_jaccard": _pair(set_jaccard),
    "index_agreement": _pair(index_agreement),
    "cell_type_nmi_difference": Comparison(("candidate", "reference", "truth"), label_nmi_difference),
    "accuracy_difference": Comparison(("candidate", "reference", "truth"), accuracy_difference),
}

# Aggregations over metrics that already exist, for criteria stated over a family.
AGGREGATIONS: dict[str, tuple[str, Callable[[list[float]], float]]] = {
    "mean_set_jaccard": ("set_jaccard", np.mean),
    "minimum_set_jaccard": ("set_jaccard", np.min),
}
