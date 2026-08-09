from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc
from _report import measure, write_report
from _shared import (
    embedding_knn_overlap,
    graph_jaccard,
    reseeded_umap_overlap,
)
from sklearn.manifold import trustworthiness
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

reference = sc.datasets.pbmc68k_reduced()
candidate = reference.copy()
rsc.get.anndata_to_GPU(candidate)

sc.pp.neighbors(reference, n_neighbors=15, n_pcs=50, random_state=0)
rsc.pp.neighbors(candidate, n_neighbors=15, n_pcs=50, algorithm="brute")
rsc.get.anndata_to_CPU(candidate)

metrics = [
    measure(
        "neighbors.distance_graph_jaccard", graph_jaccard(reference.obsp["distances"], candidate.obsp["distances"])
    ),
    measure(
        "neighbors.connectivity_graph_jaccard",
        graph_jaccard(reference.obsp["connectivities"], candidate.obsp["connectivities"]),
    ),
]

sc.tl.umap(reference, random_state=0)
rsc.tl.umap(candidate, random_state=0)

cpu_trustworthiness = trustworthiness(reference.obsm["X_pca"], reference.obsm["X_umap"], n_neighbors=15)
gpu_trustworthiness = trustworthiness(candidate.obsm["X_pca"], candidate.obsm["X_umap"], n_neighbors=15)
cross_overlap = embedding_knn_overlap(reference.obsm["X_umap"], candidate.obsm["X_umap"])

# UMAP is stochastic, so the only defensible reference for a CPU-vs-GPU overlap is
# how far the CPU reference is from itself. The baseline reruns the CPU embedding at
# other seeds; the gated criterion is that the GPU embedding is no further away than
# a reseeded CPU run is.
baseline_overlap = reseeded_umap_overlap(
    reference,
    reference.obsm["X_umap"],
    lambda adata, seed: sc.tl.umap(adata, random_state=seed),
)


# Only differences gate. Trustworthiness scores a single embedding against its own input,
# so it is recorded as evidence rather than asserted; the criterion that remains is whether
# the GPU embedding is no further from the CPU one than a reseeded CPU run is.
metrics.extend(
    [
        measure("umap.cpu.trustworthiness", cpu_trustworthiness),
        measure("umap.gpu.trustworthiness", gpu_trustworthiness),
        measure("umap.cross_embedding_knn_overlap", cross_overlap),
        measure("umap.trustworthiness_difference", abs(cpu_trustworthiness - gpu_trustworthiness)),
        measure("umap.cpu_reseeded_knn_overlap", baseline_overlap),
        measure("umap.cross_embedding_overlap_vs_cpu_baseline", cross_overlap - baseline_overlap),
    ]
)

sc.tl.leiden(reference, resolution=0.7, random_state=0, key_added="cpu_leiden", flavor="igraph")
rsc.tl.leiden(candidate, resolution=0.7, random_state=0, key_added="gpu_leiden")
metrics.extend(
    [
        measure(
            "leiden.adjusted_rand_index", adjusted_rand_score(reference.obs["cpu_leiden"], candidate.obs["gpu_leiden"])
        ),
        measure(
            "leiden.normalized_mutual_information",
            normalized_mutual_info_score(reference.obs["cpu_leiden"], candidate.obs["gpu_leiden"]),
        ),
    ]
)

write_report("scanpy_core_graphs_embeddings", "pbmc68k_reduced", "stochastic", metrics)
