from __future__ import annotations

import numpy as np
from sklearn.neighbors import NearestNeighbors


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


def ranked_names(adata, key: str, group: str, n_genes: int = 50) -> list[str]:
    names = adata.uns[key]["names"]
    if getattr(names.dtype, "names", None):
        return [str(name) for name in names[group][:n_genes]]
    group_index = list(adata.obs[adata.uns[key]["params"]["groupby"]].cat.categories).index(group)
    return [str(name) for name in names[:n_genes, group_index]]
