from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

import rapids_singlecell as rsc
import scanpy as sc
from _report import measure, write_report
from arrays import capture


def gpu_copy(adata):
    candidate = adata.copy()
    rsc.get.anndata_to_GPU(candidate)
    return candidate


METHOD = "scanpy_core_preprocessing"
metrics = []
counts = sc.datasets.pbmc3k()

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.filter_cells(cpu, min_genes=200)
rsc.pp.filter_cells(gpu, min_genes=200)
metrics.append(measure("filter_cells.obs_name_agreement", float(cpu.obs_names.equals(gpu.obs_names))))

sc.pp.filter_genes(cpu, min_cells=3)
rsc.pp.filter_genes(gpu, min_cells=3)
metrics.append(measure("filter_genes.var_name_agreement", float(cpu.var_names.equals(gpu.var_names))))

cpu = counts.copy()
gpu = gpu_copy(counts)
for adata in (cpu, gpu):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(cpu, qc_vars=["mt"], log1p=True, percent_top=False, inplace=True)
rsc.pp.calculate_qc_metrics(gpu, qc_vars=["mt"], log1p=True)
qc_columns = ["n_genes_by_counts", "total_counts", "total_counts_mt", "pct_counts_mt"]
for key in qc_columns:
    capture(METHOD, f"calculate_qc_metrics.{key}", cpu.obs[key], gpu.obs[key])

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.normalize_total(cpu, target_sum=10_000)
rsc.pp.normalize_total(gpu, target_sum=10_000)
rsc.get.anndata_to_CPU(gpu)
# The absolute error, the float32 ULP at the data's magnitude, and the scale-relative
# error are all recorded beside the criterion, because an absolute reading of this
# comparison is what made it look like a disagreement.
capture(METHOD, "normalize_total", cpu.X, gpu.X)

sc.pp.log1p(cpu)
rsc.get.anndata_to_GPU(gpu)
rsc.pp.log1p(gpu)
rsc.get.anndata_to_CPU(gpu)
capture(METHOD, "log1p", cpu.X, gpu.X)

filtered = counts.copy()
sc.pp.filter_genes(filtered, min_cells=3)
cpu = filtered.copy()
gpu = gpu_copy(filtered)
sc.experimental.pp.normalize_pearson_residuals(cpu)
rsc.pp.normalize_pearson_residuals(gpu)
rsc.get.anndata_to_CPU(gpu)
capture(METHOD, "normalize_pearson_residuals", cpu.X, gpu.X)

prepared = counts.copy()
sc.pp.calculate_qc_metrics(prepared, percent_top=None, inplace=True)
sc.pp.filter_genes(prepared, min_cells=3)
sc.pp.normalize_total(prepared, target_sum=10_000)
sc.pp.log1p(prepared)
sc.pp.highly_variable_genes(prepared, flavor="cell_ranger", n_top_genes=2_000)
prepared = prepared[:, prepared.var.highly_variable].copy()

cpu = prepared.copy()
gpu = gpu_copy(prepared)
sc.pp.scale(cpu, max_value=10)
rsc.pp.scale(gpu, max_value=10)
rsc.get.anndata_to_CPU(gpu)
capture(METHOD, "scale", cpu.X, gpu.X)

regression_input = prepared.copy()
cpu = regression_input.copy()
gpu = gpu_copy(regression_input)
sc.pp.regress_out(cpu, keys=["total_counts"])
rsc.pp.regress_out(gpu, keys=["total_counts"])
rsc.get.anndata_to_CPU(gpu)
capture(METHOD, "regress_out", cpu.X, gpu.X)

gene_list = ["CD3E", "CD8A", "IL7R", "MS4A1", "CD79A", "LYZ", "NKG7", "GNLY"]
score_input = counts.copy()
sc.pp.filter_genes(score_input, min_cells=3)
sc.pp.normalize_total(score_input, target_sum=10_000)
sc.pp.log1p(score_input)
cpu = score_input.copy()
gpu = gpu_copy(score_input)
sc.tl.score_genes(cpu, gene_list=gene_list, score_name="marker_score", random_state=0)
rsc.tl.score_genes(gpu, gene_list=gene_list, score_name="marker_score", random_state=0)
capture(METHOD, "score_genes", cpu.obs["marker_score"], gpu.obs["marker_score"])

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.sqrt(cpu)
rsc.pp.sqrt(gpu)
rsc.get.anndata_to_CPU(gpu)
capture(METHOD, "sqrt", cpu.X, gpu.X)

write_report(METHOD, "pbmc3k", "deterministic", metrics)
