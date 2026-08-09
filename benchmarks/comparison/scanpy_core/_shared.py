from __future__ import annotations

import numpy as np
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


def dense(value) -> np.ndarray:
    if sparse.issparse(value):
        return value.toarray()
    return np.asarray(value)


def max_abs(left, right) -> float:
    left_array = dense(left).astype(float, copy=False)
    right_array = dense(right).astype(float, copy=False)
    return float(np.nanmax(np.abs(left_array - right_array)))


def max_rel(left, right) -> float:
    """Largest absolute difference expressed relative to the reference scale.

    An absolute tolerance is not a meaningful criterion when the compared values
    span several orders of magnitude: in float32 the reference itself cannot be
    represented more precisely than one ULP at its largest value. Normalizing by
    that largest value gives a criterion that a correct float32 implementation
    can actually satisfy.
    """
    left_array = dense(left).astype(float, copy=False)
    right_array = dense(right).astype(float, copy=False)
    scale = float(np.nanmax(np.abs(left_array)))
    if scale == 0.0:
        return 0.0
    return float(np.nanmax(np.abs(left_array - right_array)) / scale)


def float32_ulp(left) -> float:
    """One float32 ULP at the largest magnitude present in the reference."""
    scale = float(np.nanmax(np.abs(dense(left).astype(float, copy=False))))
    return float(np.spacing(np.float32(scale)))


def pearson(left, right) -> float:
    left_array = dense(left).astype(float, copy=False).ravel()
    right_array = dense(right).astype(float, copy=False).ravel()
    mask = np.isfinite(left_array) & np.isfinite(right_array)
    if mask.sum() < 2 or np.std(left_array[mask]) == 0 or np.std(right_array[mask]) == 0:
        return float(
            left_array[mask].shape == right_array[mask].shape and np.allclose(left_array[mask], right_array[mask])
        )
    return float(np.corrcoef(left_array[mask], right_array[mask])[0, 1])


def jaccard(left, right) -> float:
    left_set = set(left)
    right_set = set(right)
    return len(left_set & right_set) / max(1, len(left_set | right_set))


def graph_jaccard(left, right) -> float:
    left = left.tocsr()
    right = right.tocsr()
    scores = []
    for index in range(left.shape[0]):
        left_neighbors = set(left.indices[left.indptr[index] : left.indptr[index + 1]]) - {index}
        right_neighbors = set(right.indices[right.indptr[index] : right.indptr[index + 1]]) - {index}
        scores.append(len(left_neighbors & right_neighbors) / max(1, len(left_neighbors | right_neighbors)))
    return float(np.mean(scores))


def embedding_knn_overlap(left, right, n_neighbors: int = 15) -> float:
    left_neighbors = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(left).kneighbors(return_distance=False)[:, 1:]
    right_neighbors = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(right).kneighbors(return_distance=False)[:, 1:]
    return float(
        np.mean(
            [
                len(set(left_row).intersection(right_row)) / n_neighbors
                for left_row, right_row in zip(left_neighbors, right_neighbors, strict=True)
            ]
        )
    )


# Seeds used to measure how far the CPU reference is from itself. UMAP is
# stochastic, so a CPU-vs-GPU overlap is only interpretable next to the overlap
# between two CPU runs that differ solely in random seed. This is recorded as
# evidence; it does not gate anything.
BASELINE_SEEDS = (1, 2)


def reseeded_umap_overlap(adata, reference_embedding, umap, *, n_neighbors: int = 15) -> float:
    """Overlap between the reference embedding and CPU reruns at other seeds.

    `umap` is called as `umap(copy, seed)` and must populate `obsm["X_umap"]`.
    The minimum over `BASELINE_SEEDS` is returned, giving the weakest
    self-consistency the reference implementation displays.
    """
    overlaps = []
    for seed in BASELINE_SEEDS:
        candidate = adata.copy()
        umap(candidate, seed)
        overlaps.append(embedding_knn_overlap(reference_embedding, candidate.obsm["X_umap"], n_neighbors=n_neighbors))
    return float(min(overlaps))


def component_abs_correlations(left, right) -> np.ndarray:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    n_components = min(left.shape[1], right.shape[1])
    return np.asarray([abs(pearson(left[:, index], right[:, index])) for index in range(n_components)])


def ranked_names(adata, key: str, group: str, n_genes: int = 50) -> list[str]:
    names = adata.uns[key]["names"]
    if getattr(names.dtype, "names", None):
        return [str(name) for name in names[group][:n_genes]]
    group_index = list(adata.obs[adata.uns[key]["params"]["groupby"]].cat.categories).index(group)
    return [str(name) for name in names[:n_genes, group_index]]


ALLCLOSE_BASIS = (
    "Manuscript Methods: deterministic operations are validated with numpy.allclose at "
    "default parameters (rtol=1e-5, atol=1e-8)."
)


def allclose_excess(left, right, *, rtol: float = 1e-5, atol: float = 1e-8) -> float:
    """Worst elementwise violation of `numpy.allclose`, as a fraction of its own envelope.

    The published validation standard for deterministic operations is `numpy.allclose` at
    its default parameters, which is a *relative* criterion: |a - b| <= atol + rtol * |b|.
    Dividing the difference by that envelope gives one scale-free number: <= 1 means the
    two arrays are `allclose`, and the value says how far past the envelope the worst
    element sits. An absolute tolerance cannot express this, because the same disagreement
    is negligible at 1e3 and fatal at 1e-3.
    """
    left_array = dense(left).astype(float, copy=False)
    right_array = dense(right).astype(float, copy=False)
    return float(np.nanmax(np.abs(left_array - right_array) / (atol + rtol * np.abs(right_array))))


def allclose_diagnosis(left, right, *, rtol: float = 1e-5, atol: float = 1e-8) -> str:
    """Explain a `numpy.allclose` verdict by locating where its envelope is spent.

    The criterion has two terms, |a - b| <= atol + rtol * |b|, and which one binds depends
    entirely on the magnitude of the element that fails. For a quantity that crosses zero
    the absolute floor binds, and at float32 precision that floor is below the noise of any
    implementation. Reporting the reference magnitude at the worst element says which of
    the two terms the verdict actually rests on.
    """
    left_array = dense(left).astype(float, copy=False)
    right_array = dense(right).astype(float, copy=False)
    envelope = atol + rtol * np.abs(right_array)
    ratio = (np.abs(left_array - right_array) / envelope).ravel()
    magnitude = float(np.abs(right_array.ravel()[int(np.nanargmax(ratio))]))
    if atol >= rtol * magnitude:
        return (
            f"The worst element sits at a reference magnitude of {magnitude:.3g}, where "
            f"numpy.allclose's absolute floor (atol={atol:g}) decides the criterion instead of its "
            f"relative term. This quantity crosses zero, and float32 cannot resolve a difference "
            f"below that floor there, so no float32 implementation can meet the default parameters "
            f"on it. The paired correlation on the same arrays says how the values agree overall."
        )
    return (
        f"The worst element sits at a reference magnitude of {magnitude:.3g}, where the relative "
        f"term (rtol={rtol:g}) decides the criterion, so this is a genuine relative disagreement."
    )
