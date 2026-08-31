from __future__ import annotations

import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
from _report import capture, write_report

METHOD = "rank_genes_groups"
adata = sc.datasets.pbmc68k_reduced()
groupby = "bulk_labels"

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
        # Top-50 gene names, for the overlap criterion.
        capture(
            METHOD,
            f"{method}.{group}.markers",
            reference=ref_names[group].iloc[:50].to_numpy(),
            candidate=gpu_names[group].iloc[:50].to_numpy(),
        )
        # Scores aligned on gene name, because the two sides rank genes in different orders
        # and correlating them positionally would compare unrelated genes.
        merged = pd.DataFrame({"gene": ref_names[group], "ref_score": ref_scores[group]}).merge(
            pd.DataFrame({"gene": gpu_names[group], "gpu_score": gpu_scores[group]}), on="gene"
        )
        capture(
            METHOD,
            f"{method}.{group}.score",
            reference=merged.ref_score.to_numpy(),
            candidate=merged.gpu_score.to_numpy(),
        )

write_report(METHOD, "scanpy.datasets.pbmc68k_reduced", "near-deterministic", [], shape=adata.shape)
