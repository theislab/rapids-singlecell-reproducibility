import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.filter_genes(adata_sc, min_cells=3)
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.experimental.pp.normalize_pearson_residuals(adata_sc)
rsc.pp.normalize_pearson_residuals(adata_rsc)
rsc.get.anndata_to_CPU(adata_rsc)
# numpy.allclose defaults, the standard the manuscript's Methods declare. Passed
# explicitly because np.testing.assert_allclose defaults to rtol=1e-7, atol=0 instead.
np.testing.assert_allclose(adata_sc.X, adata_rsc.X, rtol=1e-5, atol=1e-8)
