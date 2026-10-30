from __future__ import annotations

import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.neighbors import KNeighborsClassifier

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report
from benchmarks.comparison.scanpy_core._shared import umap_baseline

METHOD = "biological_pipeline_pbmc3k"


def pipeline(source, lib, *, leiden: dict):
    adata = source.copy()
    if lib is rsc:
        rsc.get.anndata_to_GPU(adata)
    lib.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=1_838)
    adata = adata[:, adata.var.highly_variable].copy()
    adata.layers["lognorm"] = adata.X.copy()
    lib.pp.scale(adata, max_value=10)
    lib.pp.pca(adata, n_comps=30, random_state=0)
    lib.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)
    lib.tl.umap(adata, random_state=0)
    lib.tl.leiden(adata, resolution=0.7, random_state=0, key_added="clusters", **leiden)
    lib.tl.rank_genes_groups(
        adata, groupby="cell_type", method="wilcoxon", n_genes=100, layer="lognorm", use_raw=False, key_added="markers"
    )
    if lib is rsc:
        rsc.get.anndata_to_CPU(adata, convert_all=True)
    return adata


processed = sc.datasets.pbmc3k_processed()
source = processed.raw.to_adata()
source.obs["cell_type"] = processed.obs["louvain"].astype("category")
cpu = pipeline(source, sc, leiden={"flavor": "igraph"})
gpu = pipeline(source, rsc, leiden={})

labels = np.asarray(source.obs["cell_type"])
train, test = next(StratifiedShuffleSplit(n_splits=1, test_size=0.35, random_state=0).split(labels, labels))


def annotate(embedding):
    classifier = KNeighborsClassifier(n_neighbors=7, weights="distance").fit(embedding[train], labels[train])
    return classifier.predict(embedding[test])


for group in cpu.obs["cell_type"].cat.categories:
    capture(
        METHOD,
        f"markers.{group}",
        reference=pd.DataFrame(cpu.uns["markers"]["names"])[group][:50],
        candidate=pd.DataFrame(gpu.uns["markers"]["names"])[group][:50],
    )
capture(
    METHOD, "highly_variable_genes.selection", reference=cpu.var_names.to_numpy(), candidate=gpu.var_names.to_numpy()
)
capture(METHOD, "pca", reference=cpu.obsm["X_pca"], candidate=gpu.obsm["X_pca"])
capture(METHOD, "umap", reference=cpu.obsm["X_umap"], candidate=gpu.obsm["X_umap"])
capture(
    METHOD,
    "clustering",
    reference=cpu.obs["clusters"],
    candidate=gpu.obs["clusters"],
    truth=labels,
)
capture(
    METHOD,
    "annotation",
    reference=annotate(cpu.obsm["X_pca"]),
    candidate=annotate(gpu.obsm["X_pca"]),
    truth=labels[test],
)

write_report(
    METHOD,
    "pbmc3k_processed raw log-expression with published cell-type labels",
    umap_baseline(cpu, gpu.obsm["X_umap"]),
    shape=source.shape,
)
