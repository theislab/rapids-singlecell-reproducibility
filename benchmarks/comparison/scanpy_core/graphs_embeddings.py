from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc
from _report import lower, write_report
from _shared import embedding_knn_overlap, graph_jaccard
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
metrics.extend(
    [
        lower(
            "umap.cpu.trustworthiness",
            trustworthiness(reference.obsm["X_pca"], reference.obsm["X_umap"], n_neighbors=15),
            0.9,
        ),
        lower(
            "umap.gpu.trustworthiness",
            trustworthiness(candidate.obsm["X_pca"], candidate.obsm["X_umap"], n_neighbors=15),
            0.9,
        ),
        lower(
            "umap.cross_embedding_knn_overlap",
            embedding_knn_overlap(reference.obsm["X_umap"], candidate.obsm["X_umap"]),
            0.65,
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
