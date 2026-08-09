from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def grouped_adata():
    rng = np.random.default_rng(42)
    blocks = [rng.normal(loc=i * 0.5, size=(30, 12)) for i in range(4)]
    x = np.vstack(blocks).astype(np.float64)
    obs = pd.DataFrame({"group": pd.Categorical([f"g{i}" for i in range(4) for _ in range(30)])})
    adata = ad.AnnData(x, obs=obs)
    adata.obsm["X_pca"] = x.copy()
    return adata


def guide_adata():
    rng = np.random.default_rng(42)
    x = np.zeros((200, 8), dtype=np.float32)
    for guide in range(8):
        x[:, guide] = np.concatenate([rng.poisson(2, 140), rng.poisson(50, 60)])
        rng.shuffle(x[:, guide])
    return ad.AnnData(
        x,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(200)]),
        var=pd.DataFrame(index=[f"guide_{i}" for i in range(8)]),
    )


def screen_adata():
    rng = np.random.default_rng(1)
    x = np.hstack(
        [
            np.vstack([np.clip(rng.normal(0, 1, (10, 10)), 0, None) for _ in range(3)]),
            np.vstack(
                [
                    np.clip(rng.normal(4, 0.7, (10, 10)), 0, None),
                    np.clip(rng.normal(4, 0.7, (10, 10)), 0, None),
                    np.clip(rng.normal(7, 0.9, (10, 10)), 0, None),
                ]
            ),
        ]
    ).astype(np.float32)
    obs = pd.DataFrame({"gene_target": ["NT"] * 10 + ["target_gene_a"] * 20}, index=[str(i) for i in range(30)])
    var = pd.DataFrame(index=[f"gene{i}" for i in range(20)])
    return ad.AnnData(sparse.csr_matrix(x), obs=obs, var=var)


def pearson(a, b):
    a, b = np.asarray(a, dtype=float).ravel(), np.asarray(b, dtype=float).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def max_abs(a, b):
    return float(np.nanmax(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


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
