import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.filter_genes(adata_sc, min_cells=3)

adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.experimental.pp.highly_variable_genes(adata_sc, flavor="pearson_residuals", n_top_genes=2000)
rsc.pp.highly_variable_genes(adata_rsc, flavor="pearson_residuals", n_top_genes=2000)
adata_sc = adata_sc[:, adata_sc.var.highly_variable].copy()
adata_rsc = adata_rsc[:, adata_rsc.var.highly_variable].copy()

np.testing.assert_array_equal(adata_sc.var_names, adata_rsc.var_names)
