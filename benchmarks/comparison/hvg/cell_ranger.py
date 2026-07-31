import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.filter_genes(adata_sc, min_cells=3)
sc.pp.normalize_total(adata_sc, target_sum=10000)
sc.pp.log1p(adata_sc)
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.pp.highly_variable_genes(adata_sc, flavor="cell_ranger", n_top_genes=2000)
rsc.pp.highly_variable_genes(adata_rsc, flavor="cell_ranger", n_top_genes=2000)
for column in ("means", "dispersions", "dispersions_norm"):
    np.testing.assert_allclose(adata_sc.var[column], adata_rsc.var[column], rtol=1e-5, atol=1e-6, equal_nan=True)
adata_sc = adata_sc[:, adata_sc.var.highly_variable].copy()
adata_rsc = adata_rsc[:, adata_rsc.var.highly_variable].copy()

np.testing.assert_array_equal(adata_sc.var_names, adata_rsc.var_names)
