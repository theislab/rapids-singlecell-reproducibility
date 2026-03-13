import rapids_singlecell as rsc
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

adata_sc = sc.datasets.pbmc3k()
sc.pp.filter_genes(adata_sc, min_cells=3)
sc.pp.normalize_total(adata_sc, target_sum=10000)
sc.pp.log1p(adata_sc)
sc.pp.highly_variable_genes(adata_sc, flavor="cell_ranger", n_top_genes=2000)
adata_sc = adata_sc[:, adata_sc.var.highly_variable].copy()
sc.pp.scale(adata_sc, max_value=10, zero_center=True)
sc.pp.pca(adata_sc)
sc.pp.neighbors(adata_sc)

adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)

sc.tl.leiden(adata_sc, resolution=0.7)
rsc.tl.leiden(adata_rsc, resolution=0.7)

rsc.get.anndata_to_CPU(adata_rsc)
# Leiden is stochastic; compare partitions via ARI and NMI
ari = adjusted_rand_score(adata_sc.obs["leiden"], adata_rsc.obs["leiden"])
nmi = normalized_mutual_info_score(adata_sc.obs["leiden"], adata_rsc.obs["leiden"])
assert ari > 0.9
assert nmi > 0.9
