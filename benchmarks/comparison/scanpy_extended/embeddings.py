from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
from sklearn.manifold import trustworthiness

from _report import lower, write_report
from _shared import knn_overlap, pbmc68k, pearson


adata = pbmc68k()
metrics = []

reference, candidate = adata.copy(), adata.copy()
sc.tl.tsne(reference, use_rep="X_pca", learning_rate=200, random_state=0)
rsc.tl.tsne(candidate, use_rep="X_pca", learning_rate=200)
for label, obj in (("cpu", reference), ("gpu", candidate)):
    metrics.append(lower(f"tsne.{label}.trustworthiness", trustworthiness(adata.obsm["X_pca"], obj.obsm["X_tsne"], n_neighbors=15), 0.90))
metrics.append(lower("tsne.cross_embedding_knn_overlap", knn_overlap(reference.obsm["X_tsne"], candidate.obsm["X_tsne"]), 0.70))

reference, candidate = adata.copy(), adata.copy()
sc.tl.diffmap(reference, n_comps=15)
rsc.tl.diffmap(candidate, n_comps=15)
component_corr = [abs(pearson(reference.obsm["X_diffmap"][:, i], candidate.obsm["X_diffmap"][:, i])) for i in range(1, 15)]
metrics.append(lower("diffmap.minimum_component_abs_correlation", min(component_corr), 0.95))

reference, candidate = adata.copy(), adata.copy()
sc.tl.draw_graph(reference, layout="fa", random_state=0)
rsc.tl.draw_graph(candidate, random_state=0)
metrics.append(lower("draw_graph.cross_embedding_knn_overlap", knn_overlap(reference.obsm["X_draw_graph_fa"], candidate.obsm["X_draw_graph_fa"]), 0.65))

reference, candidate = adata.copy(), adata.copy()
sc.tl.embedding_density(reference, basis="umap")
rsc.tl.embedding_density(candidate, basis="umap")
metrics.append(lower("embedding_density.pearson_correlation", pearson(reference.obs["umap_density"], candidate.obs["umap_density"]), 0.999))

write_report("embeddings_extended", "scanpy.datasets.pbmc68k_reduced", "stochastic", metrics)
