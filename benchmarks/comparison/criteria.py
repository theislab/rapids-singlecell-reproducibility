"""Every criterion the suite applies. The only place a threshold is written down.

The comparison scripts measure and record; nothing in them decides. `evaluate.py` joins
their measurements against this file, which means a threshold can be reviewed, argued
about or changed without a GPU and without re-measuring anything.

Rules are keyed by method group, because a metric name is only unique within one:
`umap.cross_embedding_knn_overlap` is `>= 0.6` on the pbmc3k pipeline and `>= 0.65` on
`scanpy_core_graphs_embeddings`, and collapsing those would silently change the science.

Each rule is `(pattern, comparison, tolerance, basis)`. `pattern` is a metric name or an
fnmatch glob for families generated in a loop. Exact names win over globs. **A measurement
with no matching rule is recorded as evidence and does not gate** — `evaluate.py` lists
those, so a new metric cannot slip in ungated unnoticed.

Thresholds are never widened to make a run green. A criterion that fails stays as it is
and carries its diagnosis.
"""

from __future__ import annotations

import fnmatch

ALLCLOSE_BASIS = (
    "Manuscript Methods: deterministic operations are validated with numpy.allclose at "
    "default parameters (rtol=1e-5, atol=1e-8)."
)

# method -> [(metric pattern, comparison, tolerance, basis)]
CRITERIA: dict[str, list[tuple[str, str, float, str]]] = {
    "bbknn_scrublet": [
        ("bbknn.connectivity_jaccard", ">=", 0.85, ""),
        ("scrublet.score_correlation", ">=", 0.95, ""),
        ("scrublet.call_agreement", ">=", 0.95, ""),
        ("scrublet_simulate_doublets.mean_library_size_relative_error", "<=", 0.02, ""),
        ("scrublet_simulate_doublets.std_library_size_relative_error", "<=", 0.05, ""),
    ],
    "biological_pipeline_pbmc3k": [
        ("highly_variable_genes.selection_jaccard", ">=", 0.98, ""),
        ("pca.minimum_component_abs_correlation", ">=", 0.95, ""),
        ("umap.cross_embedding_knn_overlap", ">=", 0.6, ""),
        ("clustering.adjusted_rand_index", ">=", 0.8, ""),
        ("clustering.normalized_mutual_information", ">=", 0.8, ""),
        ("clustering.cell_type_nmi_difference", "<=", 0.05, ""),
        ("markers.mean_top50_jaccard", ">=", 0.85, ""),
        ("markers.minimum_top50_jaccard", ">=", 0.7, ""),
        ("annotation.cpu_gpu_agreement", ">=", 0.9, ""),
        ("annotation.accuracy_difference", "<=", 0.05, ""),
    ],
    "calculate_niche": [
        ("neighborhood.adjusted_rand_index", ">=", 0.9, ""),
        ("neighborhood.normalized_mutual_information", ">=", 0.9, ""),
        ("neighborhood.cluster_count_difference", "<=", 1.0, ""),
        ("utag.adjusted_rand_index", ">=", 0.85, ""),
        ("utag.normalized_mutual_information", ">=", 0.85, ""),
        ("utag.cluster_count_difference", "<=", 1.0, ""),
        ("cellcharter.adjusted_rand_index", ">=", 0.8, ""),
        ("cellcharter.normalized_mutual_information", ">=", 0.8, ""),
        ("cellcharter.cluster_count_difference", "<=", 1.0, ""),
    ],
    "clustering_extended": [
        ("louvain.adjusted_rand_index", ">=", 0.8, ""),
        ("louvain.normalized_mutual_information", ">=", 0.8, ""),
        ("kmeans.adjusted_rand_index", ">=", 0.8, ""),
        ("kmeans.normalized_mutual_information", ">=", 0.8, ""),
    ],
    "co_occurrence": [
        ("interval.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("occurrence.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("occurrence.mean_abs_error", "<=", 1e-06, ""),
        ("occurrence.pearson_correlation", ">=", 0.99999, ""),
    ],
    "decoupler_methods": [
        ("*.allclose_excess", "<=", 1.0, ""),  # x9
        ("*.pearson_correlation", ">=", 0.999, ""),  # x9
    ],
    "distance": [
        ("*.pairwise_allclose_excess", "<=", 1.0, ""),  # x9
        ("*.pairwise_correlation", ">=", 0.9999, ""),  # x9
    ],
    "embeddings_extended": [
        ("tsne.cpu.trustworthiness", ">=", 0.9, ""),
        ("tsne.gpu.trustworthiness", ">=", 0.9, ""),
        ("tsne.cross_embedding_knn_overlap", ">=", 0.45, ""),
        ("diffmap.minimum_component_abs_correlation", ">=", 0.95, ""),
        ("draw_graph.cross_embedding_knn_overlap", ">=", 0.6, ""),
        ("embedding_density.pearson_correlation", ">=", 0.999, ""),
    ],
    "guide_assignment": [
        ("assign_by_threshold.assignment_agreement", ">=", 1.0, ""),
        ("assign_to_max_guide.assignment_agreement", ">=", 1.0, ""),
        ("assign_mixture_model.assignment_agreement", ">=", 0.9, ""),
    ],
    "ingest_cell_cycle": [
        ("ingest.label_agreement", ">=", 0.98, ""),
        ("ingest.pca_correlation", ">=", 0.999, ""),
        ("ingest.umap_knn_overlap", ">=", 0.45, ""),
        ("score_genes_cell_cycle.S_score_correlation", ">=", 0.999, ""),
        ("score_genes_cell_cycle.G2M_score_correlation", ">=", 0.999, ""),
        ("score_genes_cell_cycle.phase_agreement", ">=", 0.99, ""),
    ],
    "ligrec": [
        ("means.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("pvalues.nan_mask_agreement", ">=", 1.0, ""),
        ("pvalues.mean_abs_error", "<=", 0.05, ""),
        ("pvalues.pearson_correlation", ">=", 0.9, ""),
    ],
    "mixscale": [
        ("mixscale.score_allclose_excess", "<=", 1.0, ""),
        ("mixscale.score_correlation", ">=", 0.9999, ""),
    ],
    "mixscape": [
        ("perturbation_signature.allclose_excess", "<=", 1.0, ""),
        ("perturbation_signature.pearson_correlation", ">=", 0.999, ""),
        ("mixscape.global_class_agreement", ">=", 0.95, ""),
        ("mixscape.p_ko_correlation", ">=", 0.95, ""),
        ("mixscape.lda_abs_correlation", ">=", 0.95, ""),
    ],
    "rank_genes_groups": [
        ("*.score_correlation", ">=", 0.98, ""),  # x30
        ("*.top50_jaccard", ">=", 0.85, ""),  # x30
    ],
    "scanpy_core_graphs_embeddings": [
        ("neighbors.distance_graph_jaccard", ">=", 0.9, ""),
        ("neighbors.connectivity_graph_jaccard", ">=", 0.85, ""),
        ("umap.cross_embedding_knn_overlap", ">=", 0.65, ""),
        ("leiden.adjusted_rand_index", ">=", 0.9, ""),
        ("leiden.normalized_mutual_information", ">=", 0.9, ""),
    ],
    "scanpy_core_harmony": [
        (
            "harmony.minimum_component_abs_correlation",
            ">=",
            0.95,
            "Manuscript Methods: Harmony maintains a Pearson correlation of >95% for all corrected principal components.",
        ),
        ("harmony.mean_component_abs_correlation", ">=", 0.98, ""),
        ("harmony.standard_deviation_max_abs_error", "<=", 0.1, ""),
    ],
    "scanpy_core_hvg_pca": [
        ("*.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),  # x13
        ("highly_variable_genes.*.selection_jaccard", ">=", 0.99, ""),  # x4
        ("pca.scores.minimum_component_abs_correlation", ">=", 0.999, ""),
        ("pca.loadings.minimum_component_abs_correlation", ">=", 0.999, ""),
    ],
    "scanpy_core_preprocessing": [
        ("filter_cells.obs_name_agreement", ">=", 1.0, ""),
        ("filter_genes.var_name_agreement", ">=", 1.0, ""),
        ("calculate_qc_metrics.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("normalize_total.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("normalize_total.pearson_correlation", ">=", 0.999999, ""),
        ("log1p.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("log1p.pearson_correlation", ">=", 0.999999, ""),
        ("normalize_pearson_residuals.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("normalize_pearson_residuals.pearson_correlation", ">=", 0.99999, ""),
        ("scale.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("scale.pearson_correlation", ">=", 0.99999, ""),
        ("regress_out.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("regress_out.pearson_correlation", ">=", 0.9999, ""),
        ("score_genes.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("score_genes.pearson_correlation", ">=", 0.9999, ""),
        ("sqrt.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
    ],
    "spatial_autocorr": [
        ("moran.I.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("moran.I.pearson_correlation", ">=", 0.999999, ""),
        ("geary.C.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("geary.C.pearson_correlation", ">=", 0.999999, ""),
    ],
    "sqrt": [
        ("X.allclose_excess", "<=", 1.0, ""),
    ],
}


def criterion_for(method: str, metric: str) -> tuple[str, float, str] | None:
    """Resolve one measurement to its criterion, or None if it only records evidence."""
    rules = CRITERIA.get(method, [])
    for pattern, comparison, tolerance, basis in rules:
        if pattern == metric:
            return comparison, tolerance, basis
    for pattern, comparison, tolerance, basis in rules:
        if fnmatch.fnmatchcase(metric, pattern):
            return comparison, tolerance, basis
    return None
