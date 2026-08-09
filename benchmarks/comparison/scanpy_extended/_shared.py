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
    return float(np.mean([len(set(x).intersection(y)) / n_neighbors for x, y in zip(ka, kb, strict=True)]))


def graph_overlap(a, b):
    a, b = a.tocsr(), b.tocsr()
    scores = []
    for i in range(a.shape[0]):
        ai = set(a.indices[a.indptr[i] : a.indptr[i + 1]]) - {i}
        bi = set(b.indices[b.indptr[i] : b.indptr[i + 1]]) - {i}
        scores.append(len(ai & bi) / max(1, len(ai | bi)))
    return float(np.mean(scores))


def allclose_excess(candidate, reference, *, rtol=1e-5, atol=1e-8):
    """Worst elementwise violation of `numpy.allclose`, as a fraction of its own envelope.

    The published validation standard for deterministic operations is `numpy.allclose` at
    its default parameters, which is a *relative* criterion: |a - b| <= atol + rtol * |b|.
    Dividing the difference by that envelope gives one scale-free number: <= 1 means the
    two arrays are `allclose`, and the value says how far past the envelope the worst
    element sits. An absolute tolerance cannot express this, because the same disagreement
    is negligible at 1e3 and fatal at 1e-3.
    """
    candidate = np.asarray(candidate, dtype=float)
    reference = np.asarray(reference, dtype=float)
    return float(np.max(np.abs(candidate - reference) / (atol + rtol * np.abs(reference))))
