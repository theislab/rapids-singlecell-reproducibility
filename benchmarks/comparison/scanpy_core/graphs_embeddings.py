from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report
from benchmarks.comparison.scanpy_core._shared import gpu_copy, umap_baseline

METHOD = "scanpy_core_graphs_embeddings"

reference = sc.datasets.pbmc68k_reduced()
candidate = gpu_copy(reference)

sc.pp.neighbors(reference, n_neighbors=15, n_pcs=50, random_state=0)
rsc.pp.neighbors(candidate, n_neighbors=15, n_pcs=50, algorithm="brute")
capture(METHOD, "neighbors.distance", reference=reference.obsp["distances"], candidate=candidate.obsp["distances"])
capture(
    METHOD,
    "neighbors.connectivity",
    reference=reference.obsp["connectivities"],
    candidate=candidate.obsp["connectivities"],
)

sc.tl.umap(reference, random_state=0)
rsc.tl.umap(candidate, random_state=0)
capture(METHOD, "umap", reference=reference.obsm["X_umap"], candidate=candidate.obsm["X_umap"])

sc.tl.leiden(reference, resolution=0.7, random_state=0, key_added="cpu_leiden", flavor="igraph")
rsc.tl.leiden(candidate, resolution=0.7, random_state=0, key_added="gpu_leiden")
capture(METHOD, "leiden", reference=reference.obs["cpu_leiden"], candidate=candidate.obs["gpu_leiden"])

write_report(METHOD, "pbmc68k_reduced", umap_baseline(reference, candidate.obsm["X_umap"]), shape=reference.shape)
