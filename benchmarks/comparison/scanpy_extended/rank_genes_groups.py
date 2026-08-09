from __future__ import annotations

import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
from _report import measure, write_report
from _shared import pbmc68k, pearson

adata = pbmc68k()
groupby = "bulk_labels"
metrics = []

for method in ("t-test", "wilcoxon", "logreg"):
    reference, candidate = adata.copy(), adata.copy()
    sc.tl.rank_genes_groups(reference, groupby=groupby, method=method, n_genes=100, key_added="de")
    if method == "logreg":
        rsc.tl.rank_genes_groups_logreg(candidate, groupby=groupby, n_genes=100, key_added="de")
    else:
        rsc.tl.rank_genes_groups(candidate, groupby=groupby, method=method, n_genes=100, key_added="de")

    ref_names = pd.DataFrame(reference.uns["de"]["names"])
    gpu_names = pd.DataFrame(candidate.uns["de"]["names"])
    ref_scores = pd.DataFrame(reference.uns["de"]["scores"])
    gpu_scores = pd.DataFrame(candidate.uns["de"]["scores"])
    for group in ref_names.columns:
        top_ref = set(ref_names[group].iloc[:50])
        top_gpu = set(gpu_names[group].iloc[:50])
        metrics.append(measure(f"{method}.{group}.top50_jaccard", len(top_ref & top_gpu) / len(top_ref | top_gpu)))
        merged = pd.DataFrame({"ref_name": ref_names[group], "ref_score": ref_scores[group]}).merge(
            pd.DataFrame({"ref_name": gpu_names[group], "gpu_score": gpu_scores[group]}), on="ref_name"
        )
        metrics.append(measure(f"{method}.{group}.score_correlation", pearson(merged.ref_score, merged.gpu_score)))

write_report("rank_genes_groups", "scanpy.datasets.pbmc68k_reduced", "near-deterministic", metrics)
