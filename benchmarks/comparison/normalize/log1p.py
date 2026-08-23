import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.pp.log1p(adata_sc)
rsc.pp.log1p(adata_rsc)
rsc.get.anndata_to_CPU(adata_rsc)
np.testing.assert_allclose(adata_sc.X.toarray(), adata_rsc.X.toarray(), rtol=1e-6, atol=1e-6)
