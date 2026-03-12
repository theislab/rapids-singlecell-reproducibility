import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.calculate_qc_metrics(adata_sc, inplace=True)
sc.pp.filter_genes(adata_sc, min_cells=3)
sc.pp.normalize_total(adata_sc, target_sum=10000)
sc.pp.log1p(adata_sc)
sc.pp.highly_variable_genes(adata_sc, flavor="cell_ranger", n_top_genes=2000)
adata_sc = adata_sc[:, adata_sc.var.highly_variable].copy()
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)
sc.pp.regress_out(adata_sc, keys=["total_counts"])
rsc.pp.regress_out(adata_rsc, keys=["total_counts"])

rsc.get.anndata_to_CPU(adata_rsc)
# Due to use of parallel matrix multiplication, the results are not exactly the same as scanpy
# 64 bit would be required to get closer results
np.testing.assert_allclose(adata_sc.X, adata_rsc.X, rtol=1e-5, atol=1e-5)
