import rapids_singlecell as rsc
import scanpy as sc
from sklearn.manifold import trustworthiness

adata_sc = sc.datasets.pbmc68k_reduced()

adata_rsc = adata_sc.copy()

sc.tl.umap(adata_sc)
rsc.tl.umap(adata_rsc)

# UMAP is stochastic; verify both embeddings faithfully represent the PCA structure
rust_sc = trustworthiness(adata_sc.obsm["X_pca"], adata_sc.obsm["X_umap"], n_neighbors=15)
trust_rsc = trustworthiness(adata_sc.obsm["X_pca"], adata_rsc.obsm["X_umap"], n_neighbors=15)
assert trust_sc > 0.9
assert trust_rsc > 0.9

