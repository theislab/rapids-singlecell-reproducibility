import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)

sc.pp.filter_cells(adata_sc, min_genes=200)
rsc.pp.filter_cells(adata_rsc, min_genes=200)

np.testing.assert_array_equal(adata_sc.obs_names, adata_rsc.obs_names)
