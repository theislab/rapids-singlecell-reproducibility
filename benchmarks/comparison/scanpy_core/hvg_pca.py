from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc
from _report import capture, write_report

METHOD = "scanpy_core_hvg_pca"


def gpu_copy(adata):
    candidate = adata.copy()
    rsc.get.anndata_to_GPU(candidate)
    return candidate


counts = sc.datasets.pbmc3k()
sc.pp.filter_genes(counts, min_cells=3)

flavor_columns = {
    "seurat": ("means", "dispersions", "dispersions_norm"),
    "cell_ranger": ("means", "dispersions", "dispersions_norm"),
    "seurat_v3": ("means", "variances", "variances_norm"),
    "pearson_residuals": ("means", "variances", "residual_variances"),
}
for flavor, columns in flavor_columns.items():
    source = counts.copy()
    if flavor in {"seurat", "cell_ranger"}:
        sc.pp.normalize_total(source, target_sum=10_000)
        sc.pp.log1p(source)
    cpu = source.copy()
    gpu = gpu_copy(source)
    if flavor == "pearson_residuals":
        sc.experimental.pp.highly_variable_genes(cpu, flavor=flavor, n_top_genes=1_000)
    else:
        sc.pp.highly_variable_genes(cpu, flavor=flavor, n_top_genes=1_000)
    rsc.pp.highly_variable_genes(gpu, flavor=flavor, n_top_genes=1_000)
    capture(
        METHOD,
        f"highly_variable_genes.{flavor}.selection",
        reference=cpu.var_names[cpu.var.highly_variable].to_numpy(),
        candidate=gpu.var_names[gpu.var.highly_variable].to_numpy(),
    )
    for column in columns:
        capture(
            METHOD,
            f"highly_variable_genes.{flavor}.{column}",
            reference=cpu.var[column].to_numpy(),
            candidate=gpu.var[column].to_numpy(),
        )

pca_input = counts.copy()
sc.pp.normalize_total(pca_input, target_sum=10_000)
sc.pp.log1p(pca_input)
sc.pp.highly_variable_genes(pca_input, flavor="cell_ranger", n_top_genes=1_000)
pca_input = pca_input[:, pca_input.var.highly_variable].copy()
sc.pp.scale(pca_input, max_value=10)
cpu = pca_input.copy()
gpu = gpu_copy(pca_input)
sc.pp.pca(cpu, n_comps=50, random_state=0)
rsc.pp.pca(gpu, n_comps=50, random_state=0)
rsc.get.anndata_to_CPU(gpu)
gpu.obsm["X_pca"] = rsc.get.X_to_CPU(gpu.obsm["X_pca"])
capture(METHOD, "pca.scores", reference=cpu.obsm["X_pca"], candidate=gpu.obsm["X_pca"])
capture(METHOD, "pca.loadings", reference=cpu.varm["PCs"], candidate=gpu.varm["PCs"])
capture(
    METHOD,
    "pca.explained_variance_ratio",
    reference=cpu.uns["pca"]["variance_ratio"],
    candidate=gpu.uns["pca"]["variance_ratio"],
)

write_report(METHOD, "pbmc3k", "deterministic", [], shape=counts.shape)
