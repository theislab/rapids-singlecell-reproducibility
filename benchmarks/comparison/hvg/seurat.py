import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.normalize_total(adata_sc, target_sum=10000)
sc.pp.log1p(adata_sc)
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.pp.highly_variable_genes(adata_sc)
rsc.pp.highly_variable_genes(adata_rsc)
adata_sc = adata_sc[:, adata_sc.var.highly_variable].copy()
adata_rsc = adata_rsc[:, adata_rsc.var.highly_variable].copy()

np.testing.assert_array_equal(adata_sc.var_names, adata_rsc.var_names)
