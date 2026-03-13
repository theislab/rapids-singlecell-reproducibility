import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
sc.pp.filter_genes(adata_sc, min_cells=3)
sc.pp.normalize_total(adata_sc, target_sum=10000)
sc.pp.log1p(adata_sc)

adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)

gene_list = ["CD3E", "CD8A", "IL7R", "MS4A1", "CD79A", "LYZ", "NKG7", "GNLY"]
sc.tl.score_genes(adata_sc, gene_list=gene_list, score_name="score")
rsc.tl.score_genes(adata_rsc, gene_list=gene_list, score_name="score")

rsc.get.anndata_to_CPU(adata_rsc)
np.testing.assert_allclose(adata_sc.obs["score"].values, adata_rsc.obs["score"].values, rtol=1e-5, atol=1e-7)
