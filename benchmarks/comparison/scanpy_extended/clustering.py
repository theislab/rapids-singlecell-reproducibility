from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc
from _report import lower, write_report
from _shared import pbmc68k
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

adata = pbmc68k()
metrics = []

reference, candidate = adata.copy(), adata.copy()
sc.tl.louvain(
    reference,
    resolution=1.0,
    random_state=0,
    key_added="cpu_louvain",
    flavor="igraph",
    directed=False,
)
rsc.tl.louvain(candidate, resolution=1.0, key_added="gpu_louvain")
metrics.extend(
    [
        lower(
            "louvain.adjusted_rand_index",
            adjusted_rand_score(reference.obs["cpu_louvain"], candidate.obs["gpu_louvain"]),
            0.80,
        ),
        lower(
            "louvain.normalized_mutual_information",
            normalized_mutual_info_score(reference.obs["cpu_louvain"], candidate.obs["gpu_louvain"]),
            0.80,
        ),
    ]
)

cpu_labels = KMeans(n_clusters=8, n_init=10, random_state=42).fit_predict(adata.obsm["X_pca"][:, :50])
candidate = adata.copy()
rsc.tl.kmeans(candidate, n_clusters=8, n_pcs=50, n_init=10, random_state=42, key_added="gpu_kmeans")
metrics.extend(
    [
        lower("kmeans.adjusted_rand_index", adjusted_rand_score(cpu_labels, candidate.obs["gpu_kmeans"]), 0.80),
        lower(
            "kmeans.normalized_mutual_information",
            normalized_mutual_info_score(cpu_labels, candidate.obs["gpu_kmeans"]),
            0.80,
        ),
    ]
)

write_report("clustering_extended", "scanpy.datasets.pbmc68k_reduced", "stochastic", metrics)
