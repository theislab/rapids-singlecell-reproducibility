import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.pp.normalize_total(adata_sc)
rsc.pp.normalize_total(adata_rsc)
rsc.get.anndata_to_CPU(adata_rsc)
# numpy.allclose defaults, the standard the manuscript's Methods declare. Passed
# explicitly because np.testing.assert_allclose defaults to rtol=1e-7, atol=0 instead.
np.testing.assert_allclose(adata_sc.X.toarray(), adata_rsc.X.toarray(), rtol=1e-5, atol=1e-8)
