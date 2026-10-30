from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report

METHOD = "embeddings_extended"
adata = sc.datasets.pbmc68k_reduced()

reference, candidate = adata.copy(), adata.copy()
sc.tl.tsne(reference, use_rep="X_pca", learning_rate=200, random_state=0)
rsc.tl.tsne(candidate, use_rep="X_pca", learning_rate=200)
capture(
    METHOD,
    "tsne",
    reference=reference.obsm["X_tsne"],
    candidate=candidate.obsm["X_tsne"],
    basis=adata.obsm["X_pca"],
)

reference, candidate = adata.copy(), adata.copy()
sc.tl.diffmap(reference, n_comps=15)
rsc.tl.diffmap(candidate, n_comps=15)
# Component 0 is the trivial constant eigenvector.
capture(
    METHOD,
    "diffmap",
    reference=reference.obsm["X_diffmap"][:, 1:15],
    candidate=candidate.obsm["X_diffmap"][:, 1:15],
)

reference, candidate = adata.copy(), adata.copy()
sc.tl.draw_graph(reference, layout="fa", random_state=0)
rsc.tl.draw_graph(candidate, random_state=0)
capture(
    METHOD,
    "draw_graph",
    reference=reference.obsm["X_draw_graph_fa"],
    candidate=candidate.obsm["X_draw_graph_fa"],
)

reference, candidate = adata.copy(), adata.copy()
sc.tl.embedding_density(reference, basis="umap")
rsc.tl.embedding_density(candidate, basis="umap")
capture(
    METHOD,
    "embedding_density",
    reference=reference.obs["umap_density"],
    candidate=candidate.obs["umap_density"],
)

write_report(METHOD, "scanpy.datasets.pbmc68k_reduced", [], shape=adata.shape)
