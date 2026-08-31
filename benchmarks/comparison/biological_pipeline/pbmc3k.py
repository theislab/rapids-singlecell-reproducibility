from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))
sys.path.insert(0, str(Path(__file__).parents[1] / "scanpy_core"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
from _shared import embedding_knn_overlap, ranked_names, reseeded_umap_overlap
from report import capture, measure, write_report
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.neighbors import KNeighborsClassifier

METHOD = "biological_pipeline_pbmc3k"


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


def jaccard_overlap(left, right) -> float:
    """Only used for the artefact table; the criterion is computed from the stored names."""
    left_set, right_set = set(left), set(right)
    return len(left_set & right_set) / max(1, len(left_set | right_set))


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

# Marker lists per cell type: the overlap and its aggregate are computed at evaluation time.
marker_rows = []
for group in cpu.obs["cell_type"].cat.categories:
    cpu_markers = ranked_names(cpu, "markers", group, n_genes=50)
    gpu_markers = ranked_names(gpu, "markers", group, n_genes=50)
    capture(METHOD, f"markers.{group}", reference=np.asarray(cpu_markers), candidate=np.asarray(gpu_markers))
    marker_rows.append({"cell_type": group, "top50_jaccard": jaccard_overlap(cpu_markers, gpu_markers)})

capture(
    METHOD, "highly_variable_genes.selection", reference=cpu.var_names.to_numpy(), candidate=gpu.var_names.to_numpy()
)
capture(METHOD, "pca", reference=cpu.obsm["X_pca"], candidate=gpu.obsm["X_pca"])
capture(
    METHOD,
    "umap",
    reference=cpu.obsm["X_umap"],
    candidate=gpu.obsm["X_umap"],
    reference_basis=cpu.obsm["X_pca"],
    candidate_basis=gpu.obsm["X_pca"],
)
# Clusters travel with the published labels they are scored against, so the per-backend
# cell-type NMI and its CPU/GPU difference are derived rather than fixed here.
capture(
    METHOD,
    "clustering",
    reference=cpu.obs["clusters"],
    candidate=gpu.obs["clusters"],
    reference_truth=cpu.obs["cell_type"],
    candidate_truth=gpu.obs["cell_type"],
)
capture(
    METHOD,
    "annotation",
    reference=cpu_prediction,
    candidate=gpu_prediction,
    truth=truth,
)

# Reseeding reruns UMAP, so this cannot come from stored arrays.
baseline_overlap = reseeded_umap_overlap(
    cpu,
    cpu.obsm["X_umap"],
    lambda adata, seed: sc.tl.umap(adata, random_state=seed),
)
cross_overlap = embedding_knn_overlap(cpu.obsm["X_umap"], gpu.obsm["X_umap"])
metrics = [
    measure("umap.cpu_reseeded_knn_overlap", baseline_overlap),
    measure("umap.cross_embedding_overlap_vs_cpu_baseline", cross_overlap - baseline_overlap),
]

output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
artifact_dir = Path(os.environ.get("EQUIVALENCE_ARTIFACT_DIR", output_dir / "artifacts"))
artifact_dir.mkdir(parents=True, exist_ok=True)
pd.DataFrame(marker_rows).to_csv(artifact_dir / "biological_pipeline_marker_overlap.csv", index=False)
plot_embeddings(cpu, gpu, artifact_dir / "biological_pipeline_umap.png")

write_report(
    METHOD,
    "pbmc3k_processed raw log-expression with published cell-type labels",
    "biological",
    metrics,
    shape=source.shape,
)
