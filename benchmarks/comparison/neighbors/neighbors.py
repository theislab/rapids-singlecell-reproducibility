import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.filter_genes(adata_sc, min_cells=3)
sc.pp.normalize_total(adata_sc, target_sum=10000)
sc.pp.log1p(adata_sc)
sc.pp.highly_variable_genes(adata_sc, flavor="cell_ranger", n_top_genes=2000)
adata_sc = adata_sc[:, adata_sc.var.highly_variable].copy()
sc.pp.scale(adata_sc, max_value=10, zero_center=True)
sc.pp.pca(adata_sc)

adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)

sc.pp.neighbors(adata_sc)
rsc.pp.neighbors(adata_rsc)

rsc.get.anndata_to_CPU(adata_rsc)


def knn_indices(distances):
    """Extract sorted neighbor indices per cell, excluding self-connections (distance=0)."""
    lil = distances.tolil()
    return [sorted(j for j, d in zip(row_j, row_d) if d > 0) for row_j, row_d in zip(lil.rows, lil.data)]


# rsc may include self-connections (distance=0) not present in scanpy output
sc_nbrs = knn_indices(adata_sc.obsp["distances"])
rsc_nbrs = knn_indices(adata_rsc.obsp["distances"])
overlap = [
    len(set(sc_row).intersection(rsc_row)) / max(len(sc_row), len(rsc_row), 1)
    for sc_row, rsc_row in zip(sc_nbrs, rsc_nbrs)
]
assert np.mean(overlap) >= 0.95
assert np.quantile(overlap, 0.05) >= 0.80
