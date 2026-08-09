from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
from _report import measure, write_report
from _shared import knn_overlap, pearson

adata = sc.datasets.pbmc68k_reduced()
reference_data = adata[::2].copy()
query = adata[1::2].copy()
sc.pp.pca(reference_data)
sc.pp.neighbors(reference_data, random_state=0)
sc.tl.umap(reference_data, random_state=0)

cpu_query, gpu_query = query.copy(), query.copy()
sc.tl.ingest(cpu_query, reference_data, obs="bulk_labels", embedding_method=("pca", "umap"))
rsc.tl.ingest(gpu_query, reference_data, obs="bulk_labels", embedding_method=("pca", "umap"), algorithm="brute")
metrics = [
    measure(
        "ingest.label_agreement",
        np.mean(cpu_query.obs["bulk_labels"].to_numpy() == gpu_query.obs["bulk_labels"].to_numpy()),
    ),
    measure("ingest.pca_correlation", pearson(cpu_query.obsm["X_pca"], gpu_query.obsm["X_pca"])),
    measure("ingest.umap_knn_overlap", knn_overlap(cpu_query.obsm["X_umap"], gpu_query.obsm["X_umap"])),
]

genes = adata.var_names.tolist()
s_genes = genes[: min(20, len(genes) // 2)]
g2m_genes = genes[min(20, len(genes) // 2) : min(40, len(genes))]
cpu_cycle, gpu_cycle = adata.copy(), adata.copy()
sc.tl.score_genes_cell_cycle(cpu_cycle, s_genes=s_genes, g2m_genes=g2m_genes, random_state=0)
rsc.tl.score_genes_cell_cycle(gpu_cycle, s_genes=s_genes, g2m_genes=g2m_genes, random_state=0)
metrics.extend(
    [
        measure(
            "score_genes_cell_cycle.S_score_correlation", pearson(cpu_cycle.obs["S_score"], gpu_cycle.obs["S_score"])
        ),
        measure(
            "score_genes_cell_cycle.G2M_score_correlation",
            pearson(cpu_cycle.obs["G2M_score"], gpu_cycle.obs["G2M_score"]),
        ),
        measure(
            "score_genes_cell_cycle.phase_agreement",
            np.mean(cpu_cycle.obs["phase"].to_numpy() == gpu_cycle.obs["phase"].to_numpy()),
        ),
    ]
)

write_report("ingest_cell_cycle", "scanpy.datasets.pbmc68k_reduced", "near-deterministic", metrics)
