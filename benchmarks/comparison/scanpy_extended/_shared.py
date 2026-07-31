from __future__ import annotations

import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors


def pbmc68k():
    return sc.datasets.pbmc68k_reduced()


def pearson(a, b):
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def max_abs(a, b):
    return float(np.nanmax(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def knn_overlap(a, b, n_neighbors=15):
    a = np.asarray(a)
    b = np.asarray(b)
    ka = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(a).kneighbors(return_distance=False)[:, 1:]
    kb = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(b).kneighbors(return_distance=False)[:, 1:]
    return float(np.mean([len(set(x).intersection(y)) / n_neighbors for x, y in zip(ka, kb)]))


def graph_overlap(a, b):
    a, b = a.tocsr(), b.tocsr()
    scores = []
    for i in range(a.shape[0]):
        ai = set(a.indices[a.indptr[i] : a.indptr[i + 1]]) - {i}
        bi = set(b.indices[b.indptr[i] : b.indptr[i + 1]]) - {i}
        scores.append(len(ai & bi) / max(1, len(ai | bi)))
    return float(np.mean(scores))

