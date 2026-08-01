from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
from _report import lower, observed_only, upper, write_report
from _shared import float32_ulp, max_abs, max_rel, pearson


def gpu_copy(adata):
    candidate = adata.copy()
    rsc.get.anndata_to_GPU(candidate)
    return candidate


metrics = []
counts = sc.datasets.pbmc3k()

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.filter_cells(cpu, min_genes=200)
rsc.pp.filter_cells(gpu, min_genes=200)
metrics.append(lower("filter_cells.obs_name_agreement", float(cpu.obs_names.equals(gpu.obs_names)), 1.0))

sc.pp.filter_genes(cpu, min_cells=3)
rsc.pp.filter_genes(gpu, min_cells=3)
metrics.append(lower("filter_genes.var_name_agreement", float(cpu.var_names.equals(gpu.var_names)), 1.0))

cpu = counts.copy()
gpu = gpu_copy(counts)
for adata in (cpu, gpu):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(cpu, qc_vars=["mt"], log1p=True, percent_top=False, inplace=True)
rsc.pp.calculate_qc_metrics(gpu, qc_vars=["mt"], log1p=True)
qc_columns = ["n_genes_by_counts", "total_counts", "total_counts_mt", "pct_counts_mt"]
metrics.append(
    upper("calculate_qc_metrics.max_abs_error", max(max_abs(cpu.obs[key], gpu.obs[key]) for key in qc_columns), 1e-5)
)

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.normalize_total(cpu, target_sum=10_000)
rsc.pp.normalize_total(gpu, target_sum=10_000)
rsc.get.anndata_to_CPU(gpu)
# This criterion fails, and is left failing on purpose. pbmc3k normalized to
# target_sum=10000 reaches ~1751, where one float32 ULP is already 1.22e-4, so no
# float32 implementation can meet a 1e-5 absolute tolerance. That is recorded as
# evidence next to the criterion rather than used to change it; see ../THRESHOLDS.md.
metrics.extend(
    [
        upper("normalize_total.max_abs_error", max_abs(cpu.X, gpu.X), 1e-5),
        observed_only(
            "normalize_total.float32_ulp_at_max",
            float32_ulp(cpu.X),
            basis="Smallest representable float32 difference at the largest normalized count.",
        ),
        observed_only(
            "normalize_total.max_rel_error",
            max_rel(cpu.X, gpu.X),
            basis="Absolute error divided by the largest normalized count, for scale context.",
        ),
        lower("normalize_total.pearson_correlation", pearson(cpu.X, gpu.X), 0.999999),
    ]
)

sc.pp.log1p(cpu)
rsc.get.anndata_to_GPU(gpu)
rsc.pp.log1p(gpu)
rsc.get.anndata_to_CPU(gpu)
metrics.extend(
    [
        upper("log1p.max_abs_error", max_abs(cpu.X, gpu.X), 1e-6),
        lower("log1p.pearson_correlation", pearson(cpu.X, gpu.X), 0.999999),
    ]
)

filtered = counts.copy()
sc.pp.filter_genes(filtered, min_cells=3)
cpu = filtered.copy()
gpu = gpu_copy(filtered)
sc.experimental.pp.normalize_pearson_residuals(cpu)
rsc.pp.normalize_pearson_residuals(gpu)
rsc.get.anndata_to_CPU(gpu)
metrics.extend(
    [
        upper("normalize_pearson_residuals.max_abs_error", max_abs(cpu.X, gpu.X), 1e-5),
        lower("normalize_pearson_residuals.pearson_correlation", pearson(cpu.X, gpu.X), 0.99999),
    ]
)

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
metrics.extend(
    [
        upper("scale.max_abs_error", max_abs(cpu.X, gpu.X), 1e-5),
        lower("scale.pearson_correlation", pearson(cpu.X, gpu.X), 0.99999),
    ]
)

regression_input = prepared.copy()
cpu = regression_input.copy()
gpu = gpu_copy(regression_input)
sc.pp.regress_out(cpu, keys=["total_counts"])
rsc.pp.regress_out(gpu, keys=["total_counts"])
rsc.get.anndata_to_CPU(gpu)
metrics.extend(
    [
        upper("regress_out.max_abs_error", max_abs(cpu.X, gpu.X), 1e-4),
        lower("regress_out.pearson_correlation", pearson(cpu.X, gpu.X), 0.9999),
    ]
)

gene_list = ["CD3E", "CD8A", "IL7R", "MS4A1", "CD79A", "LYZ", "NKG7", "GNLY"]
score_input = counts.copy()
sc.pp.filter_genes(score_input, min_cells=3)
sc.pp.normalize_total(score_input, target_sum=10_000)
sc.pp.log1p(score_input)
cpu = score_input.copy()
gpu = gpu_copy(score_input)
sc.tl.score_genes(cpu, gene_list=gene_list, score_name="marker_score", random_state=0)
rsc.tl.score_genes(gpu, gene_list=gene_list, score_name="marker_score", random_state=0)
metrics.extend(
    [
        upper("score_genes.max_abs_error", max_abs(cpu.obs["marker_score"], gpu.obs["marker_score"]), 1e-5),
        lower("score_genes.pearson_correlation", pearson(cpu.obs["marker_score"], gpu.obs["marker_score"]), 0.9999),
    ]
)

cpu = counts.copy()
gpu = gpu_copy(counts)
sc.pp.sqrt(cpu)
rsc.pp.sqrt(gpu)
rsc.get.anndata_to_CPU(gpu)
metrics.append(upper("sqrt.max_abs_error", max_abs(cpu.X, gpu.X), 1e-6))

write_report("scanpy_core_preprocessing", "pbmc3k", "deterministic", metrics)
