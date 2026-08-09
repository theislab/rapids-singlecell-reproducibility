from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc
from _report import lower, observed_only, write_report
from _shared import (
    BASELINE_SEEDS,
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
    lower(
        "neighbors.distance_graph_jaccard",
        graph_jaccard(reference.obsp["distances"], candidate.obsp["distances"]),
        0.9,
    ),
    lower(
        "neighbors.connectivity_graph_jaccard",
        graph_jaccard(reference.obsp["connectivities"], candidate.obsp["connectivities"]),
        0.85,
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

QUALITY_BASIS = (
    "An absolute quality score of one implementation against its own input, so it does not test "
    "CPU/GPU agreement. The paired difference criterion is what gates."
)
OVERLAP_DIAGNOSIS = (
    "UMAP is stochastic, and this threshold asks for more agreement than the CPU reference shows "
    "against itself. See `umap.cpu_reseeded_knn_overlap`, the same measurement with only the seed "
    "changed, and `umap.cross_embedding_overlap_vs_cpu_baseline`."
)

# Only differences gate. Trustworthiness scores a single embedding against its own input,
# so it is recorded as evidence rather than asserted; the criterion that remains is whether
# the GPU embedding is no further from the CPU one than a reseeded CPU run is.
metrics.extend(
    [
        observed_only("umap.cpu.trustworthiness", cpu_trustworthiness, basis=QUALITY_BASIS),
        observed_only("umap.gpu.trustworthiness", gpu_trustworthiness, basis=QUALITY_BASIS),
        lower("umap.cross_embedding_knn_overlap", cross_overlap, 0.65) | {"diagnosis": OVERLAP_DIAGNOSIS},
        observed_only(
            "umap.trustworthiness_difference",
            abs(cpu_trustworthiness - gpu_trustworthiness),
            basis="CPU/GPU embedding-quality gap, for comparison against the CPU seed-to-seed spread.",
        ),
        observed_only(
            "umap.cpu_reseeded_knn_overlap",
            baseline_overlap,
            basis=f"Weakest CPU-vs-CPU overlap over seeds {BASELINE_SEEDS}; how far the reference is from itself.",
        ),
        observed_only(
            "umap.cross_embedding_overlap_vs_cpu_baseline",
            cross_overlap - baseline_overlap,
            basis="CPU-vs-GPU overlap minus the reseeded-CPU baseline; negative means worse than reseeding.",
        ),
    ]
)

sc.tl.leiden(reference, resolution=0.7, random_state=0, key_added="cpu_leiden", flavor="igraph")
rsc.tl.leiden(candidate, resolution=0.7, random_state=0, key_added="gpu_leiden")
metrics.extend(
    [
        lower(
            "leiden.adjusted_rand_index",
            adjusted_rand_score(reference.obs["cpu_leiden"], candidate.obs["gpu_leiden"]),
            0.9,
        ),
        lower(
            "leiden.normalized_mutual_information",
            normalized_mutual_info_score(reference.obs["cpu_leiden"], candidate.obs["gpu_leiden"]),
            0.9,
        ),
    ]
)

write_report("scanpy_core_graphs_embeddings", "pbmc68k_reduced", "stochastic", metrics)
