from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "scanpy_core"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
from _report import lower, upper, write_report
from _shared import component_abs_correlations, embedding_knn_overlap, jaccard, ranked_names
from sklearn.manifold import trustworthiness
from sklearn.metrics import accuracy_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.neighbors import KNeighborsClassifier


def cpu_pipeline(source):
    adata = source.copy()
    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=1_838)
    adata = adata[:, adata.var.highly_variable].copy()
    adata.layers["lognorm"] = adata.X.copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=30, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)
    sc.tl.umap(adata, random_state=0)
    sc.tl.leiden(adata, resolution=0.7, random_state=0, key_added="clusters", flavor="igraph")
    sc.tl.rank_genes_groups(
        adata,
        groupby="cell_type",
        method="wilcoxon",
        n_genes=100,
        layer="lognorm",
        use_raw=False,
        key_added="markers",
    )
    return adata


def gpu_pipeline(source):
    adata = source.copy()
    rsc.get.anndata_to_GPU(adata)
    rsc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=1_838)
    adata = adata[:, adata.var.highly_variable].copy()
    adata.layers["lognorm"] = adata.X.copy()
    rsc.pp.scale(adata, max_value=10)
    rsc.pp.pca(adata, n_comps=30, random_state=0)
    rsc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0, algorithm="brute")
    rsc.tl.umap(adata, random_state=0)
    rsc.tl.leiden(adata, resolution=0.7, random_state=0, key_added="clusters")
    rsc.tl.rank_genes_groups(
        adata,
        groupby="cell_type",
        method="wilcoxon",
        n_genes=100,
        layer="lognorm",
        use_raw=False,
        key_added="markers",
    )
    rsc.get.anndata_to_CPU(adata, convert_all=True)
    adata.obsm["X_pca"] = rsc.get.X_to_CPU(adata.obsm["X_pca"])
    adata.obsm["X_umap"] = rsc.get.X_to_CPU(adata.obsm["X_umap"])
    return adata


def annotate(embedding, labels):
    splitter = StratifiedShuffleSplit(n_splits=1, test_size=0.35, random_state=0)
    train, test = next(splitter.split(embedding, labels))
    classifier = KNeighborsClassifier(n_neighbors=7, weights="distance")
    classifier.fit(embedding[train], labels[train])
    return test, classifier.predict(embedding[test])


def plot_embeddings(cpu, gpu, output: Path) -> None:
    categories = list(cpu.obs["cell_type"].cat.categories)
    palette = dict(zip(categories, plt.get_cmap("tab10").colors, strict=False))
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for axis, adata, title in zip(axes, (cpu, gpu), ("CPU (Scanpy)", "GPU (rapids-singlecell)"), strict=True):
        for category in categories:
            mask = np.asarray(adata.obs["cell_type"] == category)
            axis.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=7,
                alpha=0.75,
                color=palette[category],
                label=category,
                rasterized=True,
            )
        axis.set_title(title)
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        axis.set_xticks([])
        axis.set_yticks([])
    handles, labels = axes[1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(figure)


processed = sc.datasets.pbmc3k_processed()
source = processed.raw.to_adata()
source.obs["cell_type"] = processed.obs["louvain"].astype("category")
cpu = cpu_pipeline(source)
gpu = gpu_pipeline(source)

cpu_test, cpu_prediction = annotate(cpu.obsm["X_pca"], np.asarray(cpu.obs["cell_type"]))
gpu_test, gpu_prediction = annotate(gpu.obsm["X_pca"], np.asarray(gpu.obs["cell_type"]))
if not np.array_equal(cpu_test, gpu_test):
    raise RuntimeError("CPU and GPU annotation splits differ despite a fixed random seed")
truth = np.asarray(cpu.obs["cell_type"])[cpu_test]

marker_overlaps = []
marker_rows = []
for group in cpu.obs["cell_type"].cat.categories:
    cpu_markers = ranked_names(cpu, "markers", group, n_genes=50)
    gpu_markers = ranked_names(gpu, "markers", group, n_genes=50)
    overlap = jaccard(cpu_markers, gpu_markers)
    marker_overlaps.append(overlap)
    marker_rows.append({"cell_type": group, "top50_jaccard": overlap})

cpu_accuracy = accuracy_score(truth, cpu_prediction)
gpu_accuracy = accuracy_score(truth, gpu_prediction)
cpu_label_nmi = normalized_mutual_info_score(cpu.obs["cell_type"], cpu.obs["clusters"])
gpu_label_nmi = normalized_mutual_info_score(gpu.obs["cell_type"], gpu.obs["clusters"])
pca_correlations = component_abs_correlations(cpu.obsm["X_pca"], gpu.obsm["X_pca"])

metrics = [
    lower(
        "highly_variable_genes.selection_jaccard",
        jaccard(cpu.var_names, gpu.var_names),
        0.98,
    ),
    lower("pca.minimum_component_abs_correlation", pca_correlations.min(), 0.95),
    lower(
        "umap.cpu.trustworthiness",
        trustworthiness(cpu.obsm["X_pca"], cpu.obsm["X_umap"], n_neighbors=15),
        0.9,
    ),
    lower(
        "umap.gpu.trustworthiness",
        trustworthiness(gpu.obsm["X_pca"], gpu.obsm["X_umap"], n_neighbors=15),
        0.9,
    ),
    lower(
        "umap.cross_embedding_knn_overlap",
        embedding_knn_overlap(cpu.obsm["X_umap"], gpu.obsm["X_umap"]),
        0.6,
    ),
    lower(
        "clustering.adjusted_rand_index",
        adjusted_rand_score(cpu.obs["clusters"], gpu.obs["clusters"]),
        0.8,
    ),
    lower(
        "clustering.normalized_mutual_information",
        normalized_mutual_info_score(cpu.obs["clusters"], gpu.obs["clusters"]),
        0.8,
    ),
    lower("clustering.cpu_cell_type_nmi", cpu_label_nmi, 0.6),
    lower("clustering.gpu_cell_type_nmi", gpu_label_nmi, 0.6),
    upper("clustering.cell_type_nmi_difference", abs(cpu_label_nmi - gpu_label_nmi), 0.05),
    lower("markers.mean_top50_jaccard", np.mean(marker_overlaps), 0.85),
    lower("markers.minimum_top50_jaccard", np.min(marker_overlaps), 0.7),
    lower("annotation.cpu_accuracy", cpu_accuracy, 0.85),
    lower("annotation.gpu_accuracy", gpu_accuracy, 0.85),
    lower("annotation.cpu_gpu_agreement", accuracy_score(cpu_prediction, gpu_prediction), 0.9),
    upper("annotation.accuracy_difference", abs(cpu_accuracy - gpu_accuracy), 0.05),
]

output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
artifact_dir = Path(os.environ.get("EQUIVALENCE_ARTIFACT_DIR", output_dir / "artifacts"))
artifact_dir.mkdir(parents=True, exist_ok=True)
pd.DataFrame(marker_rows).to_csv(artifact_dir / "biological_pipeline_marker_overlap.csv", index=False)
plot_embeddings(cpu, gpu, artifact_dir / "biological_pipeline_umap.png")

write_report(
    "biological_pipeline_pbmc3k",
    "pbmc3k_processed raw log-expression with published cell-type labels",
    "biological",
    metrics,
)
