from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc

from benchmarks.comparison.comparisons import cross_embedding_knn_overlap
from benchmarks.comparison.report import measure


def gpu_copy(adata):
    candidate = adata.copy()
    rsc.get.anndata_to_GPU(candidate)
    return candidate


def umap_baseline(cpu, gpu_umap) -> list[dict]:
    """GPU-vs-CPU UMAP overlap against the CPU's own overlap with reseeded runs."""
    reseeded = []
    for seed in (1, 2):
        rerun = cpu.copy()
        sc.tl.umap(rerun, random_state=seed)
        reseeded.append(cross_embedding_knn_overlap(rerun.obsm["X_umap"], cpu.obsm["X_umap"]))
    baseline = min(reseeded)
    return [
        measure("umap.cpu_reseeded_knn_overlap", baseline),
        measure(
            "umap.cross_embedding_overlap_vs_cpu_baseline",
            cross_embedding_knn_overlap(gpu_umap, cpu.obsm["X_umap"]) - baseline,
        ),
    ]
