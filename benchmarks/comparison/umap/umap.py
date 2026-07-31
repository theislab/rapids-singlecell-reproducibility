import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
from sklearn.manifold import trustworthiness
from sklearn.neighbors import NearestNeighbors

adata_sc = sc.datasets.pbmc68k_reduced()

adata_rsc = adata_sc.copy()

sc.tl.umap(adata_sc)
rsc.tl.umap(adata_rsc)

# UMAP is stochastic; verify both embeddings faithfully represent the PCA structure
rust_sc = trustworthiness(adata_sc.obsm["X_pca"], adata_sc.obsm["X_umap"], n_neighbors=15)
trust_rsc = trustworthiness(adata_sc.obsm["X_pca"], adata_rsc.obsm["X_umap"], n_neighbors=15)
assert trust_sc > 0.9
assert trust_rsc > 0.9


# Compare the embeddings to each other through the local neighborhoods they imply.
# This is invariant to rotation, reflection, translation, and global scale.
def embedding_neighbors(embedding, n_neighbors=15):
    return NearestNeighbors(n_neighbors=n_neighbors + 1).fit(embedding).kneighbors(return_distance=False)[:, 1:]


sc_neighbors = embedding_neighbors(adata_sc.obsm["X_umap"])
rsc_neighbors = embedding_neighbors(adata_rsc.obsm["X_umap"])
overlap = np.mean(
    [len(set(cpu).intersection(gpu)) / len(cpu) for cpu, gpu in zip(sc_neighbors, rsc_neighbors, strict=True)]
)
assert overlap >= 0.70
