from __future__ import annotations

from typing import TYPE_CHECKING

import squidpy as sq

if TYPE_CHECKING:
    import numpy as np

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
