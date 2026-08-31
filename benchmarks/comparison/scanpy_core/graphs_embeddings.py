from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc
from _report import capture, measure, write_report
from _shared import BASELINE_SEEDS, embedding_knn_overlap, reseeded_umap_overlap

METHOD = "scanpy_core_graphs_embeddings"

reference = sc.datasets.pbmc68k_reduced()
candidate = reference.copy()
rsc.get.anndata_to_GPU(candidate)

sc.pp.neighbors(reference, n_neighbors=15, n_pcs=50, random_state=0)
rsc.pp.neighbors(candidate, n_neighbors=15, n_pcs=50, algorithm="brute")
rsc.get.anndata_to_CPU(candidate)

capture(METHOD, "neighbors.distance", reference=reference.obsp["distances"], candidate=candidate.obsp["distances"])
capture(
    METHOD,
    "neighbors.connectivity",
    reference=reference.obsp["connectivities"],
    candidate=candidate.obsp["connectivities"],
)

sc.tl.umap(reference, random_state=0)
rsc.tl.umap(candidate, random_state=0)

# The embeddings and the basis each was built from are both stored, so trustworthiness —
# which scores an embedding against its own input — and the cross-embedding overlap are
# computed at evaluation time rather than fixed here.
capture(
    METHOD,
    "umap",
    reference=reference.obsm["X_umap"],
    candidate=candidate.obsm["X_umap"],
    reference_basis=reference.obsm["X_pca"],
    candidate_basis=candidate.obsm["X_pca"],
)

# UMAP is stochastic, so the only defensible reference for a CPU-vs-GPU overlap is how far
# the CPU reference is from itself. This genuinely reruns UMAP at other seeds, so it cannot
# be recovered from stored arrays and stays a measurement.
baseline_overlap = reseeded_umap_overlap(
    reference,
    reference.obsm["X_umap"],
    lambda adata, seed: sc.tl.umap(adata, random_state=seed),
)
cross_overlap = embedding_knn_overlap(reference.obsm["X_umap"], candidate.obsm["X_umap"])
metrics = [
    measure("umap.cpu_reseeded_knn_overlap", baseline_overlap),
    measure("umap.cross_embedding_overlap_vs_cpu_baseline", cross_overlap - baseline_overlap),
]

sc.tl.leiden(reference, resolution=0.7, random_state=0, key_added="cpu_leiden", flavor="igraph")
rsc.tl.leiden(candidate, resolution=0.7, random_state=0, key_added="gpu_leiden")
capture(METHOD, "leiden", reference=reference.obs["cpu_leiden"], candidate=candidate.obs["gpu_leiden"])

write_report(METHOD, "pbmc68k_reduced", "stochastic", metrics, shape=reference.shape)
