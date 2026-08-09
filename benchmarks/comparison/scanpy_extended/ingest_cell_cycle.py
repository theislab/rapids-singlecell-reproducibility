from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc
from _report import capture, write_report

METHOD = "ingest_cell_cycle"

adata = sc.datasets.pbmc68k_reduced()
reference_data = adata[::2].copy()
query = adata[1::2].copy()
sc.pp.pca(reference_data)
sc.pp.neighbors(reference_data, random_state=0)
sc.tl.umap(reference_data, random_state=0)

cpu_query, gpu_query = query.copy(), query.copy()
sc.tl.ingest(cpu_query, reference_data, obs="bulk_labels", embedding_method=("pca", "umap"))
rsc.tl.ingest(gpu_query, reference_data, obs="bulk_labels", embedding_method=("pca", "umap"), algorithm="brute")
capture(METHOD, "ingest.label", reference=cpu_query.obs["bulk_labels"], candidate=gpu_query.obs["bulk_labels"])
capture(METHOD, "ingest.pca", reference=cpu_query.obsm["X_pca"], candidate=gpu_query.obsm["X_pca"])
capture(METHOD, "ingest.umap", reference=cpu_query.obsm["X_umap"], candidate=gpu_query.obsm["X_umap"])

genes = adata.var_names.tolist()
s_genes = genes[: min(20, len(genes) // 2)]
g2m_genes = genes[min(20, len(genes) // 2) : min(40, len(genes))]
cpu_cycle, gpu_cycle = adata.copy(), adata.copy()
sc.tl.score_genes_cell_cycle(cpu_cycle, s_genes=s_genes, g2m_genes=g2m_genes, random_state=0)
rsc.tl.score_genes_cell_cycle(gpu_cycle, s_genes=s_genes, g2m_genes=g2m_genes, random_state=0)
capture(METHOD, "score_genes_cell_cycle.S", reference=cpu_cycle.obs["S_score"], candidate=gpu_cycle.obs["S_score"])
capture(
    METHOD, "score_genes_cell_cycle.G2M", reference=cpu_cycle.obs["G2M_score"], candidate=gpu_cycle.obs["G2M_score"]
)
capture(METHOD, "score_genes_cell_cycle.phase", reference=cpu_cycle.obs["phase"], candidate=gpu_cycle.obs["phase"])

write_report(METHOD, "scanpy.datasets.pbmc68k_reduced", "near-deterministic", [])
