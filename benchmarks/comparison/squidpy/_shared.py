from __future__ import annotations

import numpy as np
import squidpy as sq

IMC_CLUSTER_KEY = "cell type"


def load_imc():
    """Load the labeled IMC example and construct one shared spatial graph."""
    adata = sq.datasets.imc()
    adata.obs[IMC_CLUSTER_KEY] = adata.obs[IMC_CLUSTER_KEY].astype("category")
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    return adata


def dataframe_values(frame) -> np.ndarray:
    """Return dense values for either a dense or pandas-sparse DataFrame."""
    try:
        frame = frame.sparse.to_dense()
    except AttributeError:
        pass
    return frame.to_numpy()



def mean_abs_error(reference, candidate) -> float:
    reference = np.asarray(reference, dtype=np.float64)
    candidate = np.asarray(candidate, dtype=np.float64)
    return float(np.nanmean(np.abs(reference - candidate)))


def pearson_correlation(reference, candidate) -> float:
    reference = np.asarray(reference, dtype=np.float64).ravel()
    candidate = np.asarray(candidate, dtype=np.float64).ravel()
    mask = np.isfinite(reference) & np.isfinite(candidate)
    return float(np.corrcoef(reference[mask], candidate[mask])[0, 1])


def allclose_excess(candidate, reference, *, rtol: float = 1e-5, atol: float = 1e-8) -> float:
    """Worst elementwise violation of `numpy.allclose`, as a fraction of its own envelope.

    The published validation standard for deterministic operations is `numpy.allclose` at
    its default parameters, which is a *relative* criterion: |a - b| <= atol + rtol * |b|.
    Dividing the difference by that envelope gives one scale-free number: <= 1 means the
    two arrays are `allclose`, and the value says how far past the envelope the worst
    element sits. An absolute tolerance cannot express this, because the same disagreement
    is negligible at 1e3 and fatal at 1e-3.
    """
    reference = np.asarray(reference, dtype=np.float64)
    candidate = np.asarray(candidate, dtype=np.float64)
    return float(np.max(np.abs(candidate - reference) / (atol + rtol * np.abs(reference))))
