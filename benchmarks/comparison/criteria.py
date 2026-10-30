"""Every criterion the suite applies, keyed by method group."""

from __future__ import annotations

import fnmatch

# method -> [(metric pattern, comparison, tolerance)]
CRITERIA: dict[str, list[tuple[str, str, float]]] = {
    "bbknn_scrublet": [
        ("bbknn.graph_jaccard", ">=", 0.85),
        ("scrublet.score.pearson_correlation", ">=", 0.95),
        ("scrublet.call.exact_agreement", ">=", 0.95),
        ("scrublet_simulate_doublets.library_size.relative_error_of_mean", "<=", 0.02),
        ("scrublet_simulate_doublets.library_size.relative_error_of_std", "<=", 0.05),
    ],
    "biological_pipeline_pbmc3k": [
        ("highly_variable_genes.selection.set_jaccard", ">=", 0.98),
        ("pca.minimum_component_abs_correlation", ">=", 0.95),
        # UMAP is stochastic: GPU-vs-CPU overlap is held to the CPU's own overlap with a reseeded run.
        ("umap.cross_embedding_overlap_vs_cpu_baseline", ">=", -0.05),
        ("clustering.adjusted_rand_index", ">=", 0.8),
        ("clustering.normalized_mutual_information", ">=", 0.8),
        ("clustering.cell_type_nmi_difference", "<=", 0.05),
        ("markers.mean_set_jaccard", ">=", 0.85),
        ("markers.minimum_set_jaccard", ">=", 0.7),
        ("annotation.exact_agreement", ">=", 0.9),
        ("annotation.accuracy_difference", "<=", 0.05),
    ],
    "co_occurrence": [
        ("occurrence.mean_abs_error", "<=", 1e-06),
        ("occurrence.pearson_correlation", ">=", 0.99999),
    ],
    "decoupler_methods": [
        ("*.pearson_correlation", ">=", 0.999),
    ],
    "distance": [
        ("*.pairwise.pearson_correlation", ">=", 0.9999),
        ("*.onesided.pearson_correlation", ">=", 0.9999),
    ],
    "embeddings_extended": [
        ("tsne.cpu.trustworthiness", ">=", 0.9),
        ("tsne.gpu.trustworthiness", ">=", 0.9),
        ("tsne.cross_embedding_knn_overlap", ">=", 0.45),
        ("diffmap.minimum_component_abs_correlation", ">=", 0.95),
        ("draw_graph.cross_embedding_knn_overlap", ">=", 0.6),
        ("embedding_density.pearson_correlation", ">=", 0.999),
    ],
    "guide_assignment": [
        ("assign_by_threshold.exact_agreement", ">=", 1.0),
        ("assign_to_max_guide.exact_agreement", ">=", 1.0),
        ("assign_mixture_model.exact_agreement", ">=", 0.9),
    ],
    "ingest_cell_cycle": [
        ("ingest.label.exact_agreement", ">=", 0.98),
        ("ingest.pca.pearson_correlation", ">=", 0.999),
        ("ingest.umap.cross_embedding_knn_overlap", ">=", 0.45),
        ("score_genes_cell_cycle.S.pearson_correlation", ">=", 0.999),
        ("score_genes_cell_cycle.G2M.pearson_correlation", ">=", 0.999),
        ("score_genes_cell_cycle.phase.exact_agreement", ">=", 0.99),
    ],
    "ligrec": [
        ("means.float32_allclose_excess", "<=", 1.0),
        ("pvalues.nan_mask_agreement", ">=", 1.0),
        ("pvalues.mean_abs_error", "<=", 0.05),
        ("pvalues.pearson_correlation", ">=", 0.9),
    ],
    "mixscale": [
        ("mixscale.score.pearson_correlation", ">=", 0.9999),
    ],
    "mixscape": [
        ("perturbation_signature.pearson_correlation", ">=", 0.999),
        ("mixscape.global_class.exact_agreement", ">=", 0.95),
        ("mixscape.p_ko.pearson_correlation", ">=", 0.95),
        ("mixscape.lda.lda_abs_correlation", ">=", 0.95),
    ],
    "rank_genes_groups": [
        ("*.score.pearson_correlation", ">=", 0.98),
        ("*.markers.set_jaccard", ">=", 0.85),
    ],
    "scanpy_core_graphs_embeddings": [
        ("neighbors.distance.graph_jaccard", ">=", 0.9),
        ("neighbors.connectivity.graph_jaccard", ">=", 0.85),
        ("umap.cross_embedding_overlap_vs_cpu_baseline", ">=", -0.05),
        ("leiden.adjusted_rand_index", ">=", 0.9),
        ("leiden.normalized_mutual_information", ">=", 0.9),
    ],
    "scanpy_core_harmony": [
        ("harmony.minimum_component_abs_correlation", ">", 0.95),
        ("harmony.mean_component_abs_correlation", ">=", 0.98),
        ("harmony.standard_deviation_max_abs_error", "<=", 0.1),
        ("harmony.relative_l2_max", "<=", 0.1),
    ],
    "scanpy_core_hvg_pca": [
        ("highly_variable_genes.seurat.selection.set_jaccard", ">=", 0.99),
        ("highly_variable_genes.seurat.means.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.seurat.dispersions.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.seurat.dispersions_norm.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.cell_ranger.selection.set_jaccard", ">=", 0.99),
        ("highly_variable_genes.cell_ranger.means.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.cell_ranger.dispersions.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.cell_ranger.dispersions_norm.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.seurat_v3.selection.set_jaccard", ">=", 0.99),
        ("highly_variable_genes.seurat_v3.means.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.seurat_v3.variances.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.seurat_v3.variances_norm.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.pearson_residuals.selection.set_jaccard", ">=", 0.99),
        ("highly_variable_genes.pearson_residuals.means.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.pearson_residuals.variances.float32_allclose_excess", "<=", 1.0),
        ("highly_variable_genes.pearson_residuals.residual_variances.float32_allclose_excess", "<=", 1.0),
        ("pca.scores.minimum_component_abs_correlation", ">=", 0.95),
        ("pca.loadings.minimum_component_abs_correlation", ">=", 0.95),
        ("pca.explained_variance_ratio.float32_allclose_excess", "<=", 1.0),
        # Magnitudes, because a principal component's sign is arbitrary.
        ("pca.scores.abs_float32_allclose_excess", "<=", 1.0),
        ("pca.loadings.abs_float32_allclose_excess", "<=", 1.0),
    ],
    "scanpy_core_preprocessing": [
        ("filter_cells.index_agreement", ">=", 1.0),
        ("filter_genes.index_agreement", ">=", 1.0),
        ("calculate_qc_metrics.*.float32_allclose_excess", "<=", 1.0),
        ("normalize_total.float32_allclose_excess", "<=", 1.0),
        ("normalize_total.pearson_correlation", ">=", 0.999999),
        ("log1p.float32_allclose_excess", "<=", 1.0),
        ("log1p.pearson_correlation", ">=", 0.999999),
        ("normalize_pearson_residuals.float32_allclose_excess", "<=", 1.0),
        ("normalize_pearson_residuals.pearson_correlation", ">=", 0.99999),
        ("scale.float32_allclose_excess", "<=", 1.0),
        ("scale.pearson_correlation", ">=", 0.99999),
        ("regress_out.pearson_correlation", ">=", 0.9999),
        ("score_genes.pearson_correlation", ">=", 0.9999),
        ("sqrt.float32_allclose_excess", "<=", 1.0),
    ],
    "spatial_autocorr": [
        ("moran.I.pearson_correlation", ">=", 0.999999),
        ("geary.C.pearson_correlation", ">=", 0.999999),
    ],
}

# Computed and reported, not gated.
EVIDENCE: dict[str, list[str]] = {
    "co_occurrence": ["occurrence.allclose_excess"],
    "decoupler_methods": ["*.allclose_excess"],
    "distance": ["*.pairwise.allclose_excess", "*.onesided.allclose_excess"],
    "mixscale": ["mixscale.score.allclose_excess"],
    "mixscape": ["perturbation_signature.allclose_excess"],
    "biological_pipeline_pbmc3k": [
        "umap.cross_embedding_knn_overlap",
        "umap.cpu_reseeded_knn_overlap",
        "markers.*.set_jaccard",
    ],
    "scanpy_core_graphs_embeddings": [
        "umap.cross_embedding_knn_overlap",
        "umap.cpu_reseeded_knn_overlap",
        "neighbors.distance.graph_exact_agreement",
        "neighbors.connectivity.graph_exact_agreement",
    ],
    "scanpy_core_preprocessing": [
        "regress_out.allclose_excess",
        "score_genes.allclose_excess",
    ],
    "spatial_autocorr": ["moran.I.allclose_excess", "geary.C.allclose_excess"],
}

# Every float32-gated criterion also reports its value at numpy's default atol, without a verdict.
for _method, _rules in CRITERIA.items():
    _defaults = [metric.replace("float32_", "") for metric, *_ in _rules if "float32_allclose_excess" in metric]
    if _defaults:
        EVIDENCE.setdefault(_method, []).extend(_defaults)


def criterion_for(method: str, metric: str) -> tuple[str, float] | None:
    """Resolve one measurement to its criterion, or None if it only records evidence."""
    for pattern, comparison, tolerance in CRITERIA.get(method, []):
        if fnmatch.fnmatchcase(metric, pattern):
            return comparison, tolerance
    return None
