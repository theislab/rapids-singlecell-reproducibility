import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.filter_genes(adata_sc, min_cells=3)
sc.pp.normalize_total(adata_sc, target_sum=10000)
sc.pp.log1p(adata_sc)
sc.pp.highly_variable_genes(adata_sc, flavor="cell_ranger", n_top_genes=2000)
adata_sc = adata_sc[:, adata_sc.var.highly_variable].copy()
sc.pp.scale(adata_sc, zero_center=True)
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.pp.pca(adata_sc)
rsc.pp.pca(adata_rsc)
# PCA Loadings
np.testing.assert_allclose(np.abs(adata_sc.varm["PCs"]), np.abs(adata_rsc.varm["PCs"]), rtol=1e-6, atol=1e-6)
# PCA Transform
np.testing.assert_allclose(np.abs(adata_sc.obsm["X_pca"]), np.abs(adata_rsc.obsm["X_pca"]), rtol=1e-6, atol=1e-6)
