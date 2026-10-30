from __future__ import annotations

import squidpy as sq

IMC_CLUSTER_KEY = "cell type"


def load_imc():
    """Load the labeled IMC example and construct one shared spatial graph."""
    adata = sq.datasets.imc()
    adata.obs[IMC_CLUSTER_KEY] = adata.obs[IMC_CLUSTER_KEY].astype("category")
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    return adata
