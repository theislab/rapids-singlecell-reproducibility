import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)

sc.pp.filter_genes(adata_sc, min_cells=3)
rsc.pp.filter_genes(adata_rsc, min_cells=3)

np.testing.assert_array_equal(adata_sc.var_names, adata_rsc.var_names)
