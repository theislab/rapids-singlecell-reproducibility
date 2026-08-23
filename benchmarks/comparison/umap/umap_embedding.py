# Named umap_embedding.py, not umap.py: a script run as `python umap/umap.py` puts its own
# directory first on sys.path, so the file would shadow the installed `umap` package and
# umap-learn's own `import umap.umap_` would fail. The sibling directory name is harmless --
# a directory without __init__.py is only a namespace portion and loses to a real package.
import rapids_singlecell as rsc
import scanpy as sc
from sklearn.manifold import trustworthiness

adata_sc = sc.datasets.pbmc68k_reduced()

adata_rsc = adata_sc.copy()

sc.tl.umap(adata_sc)
rsc.tl.umap(adata_rsc)

# UMAP is stochastic; verify both embeddings faithfully represent the PCA structure
trust_sc = trustworthiness(adata_sc.obsm["X_pca"], adata_sc.obsm["X_umap"], n_neighbors=15)
trust_rsc = trustworthiness(adata_sc.obsm["X_pca"], adata_rsc.obsm["X_umap"], n_neighbors=15)
assert trust_sc > 0.9
assert trust_rsc > 0.9

