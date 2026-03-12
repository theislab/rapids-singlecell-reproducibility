import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.filter_genes(adata_sc, min_cells=3)
sc.pp.normalize_total(adata_sc, target_sum=10000)
sc.pp.log1p(adata_sc)
sc.pp.highly_variable_genes(adata_sc, flavor="cell_ranger", n_top_genes=2000)
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.pp.scale(adata_sc, max_value=10, zero_center=True)
rsc.pp.scale(adata_rsc, max_value=10, zero_center=True)
rsc.get.anndata_to_CPU(adata_rsc)
np.testing.assert_allclose(adata_sc.X, adata_rsc.X, rtol=1e-6, atol=1e-6)
