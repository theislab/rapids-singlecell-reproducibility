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
