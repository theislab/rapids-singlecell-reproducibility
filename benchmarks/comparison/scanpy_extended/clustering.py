from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc
from _report import capture, write_report
from sklearn.cluster import KMeans

METHOD = "clustering_extended"

adata = sc.datasets.pbmc68k_reduced()
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
capture(METHOD, "louvain", reference=reference.obs["cpu_louvain"], candidate=candidate.obs["gpu_louvain"])

cpu_labels = KMeans(n_clusters=8, n_init=10, random_state=42).fit_predict(adata.obsm["X_pca"][:, :50])
candidate = adata.copy()
rsc.tl.kmeans(candidate, n_clusters=8, n_pcs=50, n_init=10, random_state=42, key_added="gpu_kmeans")
capture(METHOD, "kmeans", reference=cpu_labels, candidate=candidate.obs["gpu_kmeans"])

write_report(METHOD, "scanpy.datasets.pbmc68k_reduced", "stochastic", metrics, shape=adata.shape)
