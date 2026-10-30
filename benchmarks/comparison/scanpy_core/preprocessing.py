from __future__ import annotations

import rapids_singlecell as rsc
import scanpy as sc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report
from benchmarks.comparison.scanpy_core._shared import gpu_copy

METHOD = "scanpy_core_preprocessing"
counts = sc.datasets.pbmc3k()

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.filter_cells(cpu, min_genes=200)
rsc.pp.filter_cells(gpu, min_genes=200)
capture(METHOD, "filter_cells", reference=cpu.obs_names.to_numpy(), candidate=gpu.obs_names.to_numpy())

sc.pp.filter_genes(cpu, min_cells=3)
rsc.pp.filter_genes(gpu, min_cells=3)
capture(METHOD, "filter_genes", reference=cpu.var_names.to_numpy(), candidate=gpu.var_names.to_numpy())

cpu = counts.copy()
gpu = gpu_copy(counts)
QC_VARS = ["mt", "ribo", "hb"]
for adata in (cpu, gpu):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(cpu, qc_vars=QC_VARS, log1p=False, percent_top=False, inplace=True)
rsc.pp.calculate_qc_metrics(gpu, qc_vars=QC_VARS, log1p=False)
qc_columns = ["n_genes_by_counts", "total_counts"]
qc_columns += [f"{stat}_{var}" for var in QC_VARS for stat in ("total_counts", "pct_counts")]
for key in qc_columns:
    capture(METHOD, f"calculate_qc_metrics.{key}", reference=cpu.obs[key], candidate=gpu.obs[key])

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.normalize_total(cpu, target_sum=10_000)
rsc.pp.normalize_total(gpu, target_sum=10_000)
capture(METHOD, "normalize_total", reference=cpu.X, candidate=gpu.X)

sc.pp.log1p(cpu)
rsc.pp.log1p(gpu)
capture(METHOD, "log1p", reference=cpu.X, candidate=gpu.X)

filtered = counts.copy()
sc.pp.filter_genes(filtered, min_cells=3)
cpu = filtered.copy()
gpu = gpu_copy(filtered)
sc.experimental.pp.normalize_pearson_residuals(cpu)
rsc.pp.normalize_pearson_residuals(gpu)
capture(METHOD, "normalize_pearson_residuals", reference=cpu.X, candidate=gpu.X)

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
capture(METHOD, "scale", reference=cpu.X, candidate=gpu.X)

cpu = prepared.copy()
gpu = gpu_copy(prepared)
sc.pp.regress_out(cpu, keys=["total_counts"])
rsc.pp.regress_out(gpu, keys=["total_counts"])
capture(METHOD, "regress_out", reference=cpu.X, candidate=gpu.X)

gene_list = ["CD3E", "CD8A", "IL7R", "MS4A1", "CD79A", "LYZ", "NKG7", "GNLY"]
score_input = counts.copy()
sc.pp.filter_genes(score_input, min_cells=3)
sc.pp.normalize_total(score_input, target_sum=10_000)
sc.pp.log1p(score_input)
cpu = score_input.copy()
gpu = gpu_copy(score_input)
sc.tl.score_genes(cpu, gene_list=gene_list, score_name="marker_score", random_state=0)
rsc.tl.score_genes(gpu, gene_list=gene_list, score_name="marker_score", random_state=0)
capture(METHOD, "score_genes", reference=cpu.obs["marker_score"], candidate=gpu.obs["marker_score"])

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.sqrt(cpu)
rsc.pp.sqrt(gpu)
capture(METHOD, "sqrt", reference=cpu.X, candidate=gpu.X)

write_report(METHOD, "pbmc3k", [], shape=counts.shape)
