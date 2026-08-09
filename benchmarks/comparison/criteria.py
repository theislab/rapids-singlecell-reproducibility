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
        ("scrublet.score.pearson_correlation", ">=", 0.95, ""),
        ("scrublet.call.exact_agreement", ">=", 0.95, ""),
        ("scrublet_simulate_doublets.library_size.relative_error_of_mean", "<=", 0.02, ""),
        ("scrublet_simulate_doublets.library_size.relative_error_of_std", "<=", 0.05, ""),
    ],
    "biological_pipeline_pbmc3k": [
        ("highly_variable_genes.selection.set_jaccard", ">=", 0.98, ""),
        ("pca.minimum_component_abs_correlation", ">=", 0.95, ""),
        ("umap.cross_embedding_knn_overlap", ">=", 0.6, ""),
        ("clustering.adjusted_rand_index", ">=", 0.8, ""),
        ("clustering.normalized_mutual_information", ">=", 0.8, ""),
        ("clustering.cell_type_nmi_difference", "<=", 0.05, ""),
        ("markers.mean_set_jaccard", ">=", 0.85, ""),
        ("markers.minimum_set_jaccard", ">=", 0.7, ""),
        ("annotation.exact_agreement", ">=", 0.9, ""),
        ("annotation.accuracy_difference", "<=", 0.05, ""),
    ],
    "calculate_niche": [
        ("neighborhood.adjusted_rand_index", ">=", 0.9, ""),
        ("neighborhood.normalized_mutual_information", ">=", 0.9, ""),
        ("*.cluster_count_difference", "<=", 1.0, ""),
        ("utag.adjusted_rand_index", ">=", 0.85, ""),
        ("utag.normalized_mutual_information", ">=", 0.85, ""),
        ("cellcharter.adjusted_rand_index", ">=", 0.8, ""),
        ("cellcharter.normalized_mutual_information", ">=", 0.8, ""),
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
        ("*.allclose_excess", "<=", 1.0, ""),
        ("*.pearson_correlation", ">=", 0.999, ""),
    ],
    "distance": [
        ("*.pairwise.allclose_excess", "<=", 1.0, ""),
        ("*.pairwise.pearson_correlation", ">=", 0.9999, ""),
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
        ("assign_by_threshold.exact_agreement", ">=", 1.0, ""),
        ("assign_to_max_guide.exact_agreement", ">=", 1.0, ""),
        ("assign_mixture_model.exact_agreement", ">=", 0.9, ""),
    ],
    "ingest_cell_cycle": [
        ("ingest.label.exact_agreement", ">=", 0.98, ""),
        ("ingest.pca.pearson_correlation", ">=", 0.999, ""),
        ("ingest.umap.cross_embedding_knn_overlap", ">=", 0.45, ""),
        ("score_genes_cell_cycle.S.pearson_correlation", ">=", 0.999, ""),
        ("score_genes_cell_cycle.G2M.pearson_correlation", ">=", 0.999, ""),
        ("score_genes_cell_cycle.phase.exact_agreement", ">=", 0.99, ""),
    ],
    "ligrec": [
        ("means.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("pvalues.nan_mask_agreement", ">=", 1.0, ""),
        ("pvalues.mean_abs_error", "<=", 0.05, ""),
        ("pvalues.pearson_correlation", ">=", 0.9, ""),
    ],
    "mixscale": [
        ("mixscale.score.allclose_excess", "<=", 1.0, ""),
        ("mixscale.score.pearson_correlation", ">=", 0.9999, ""),
    ],
    "mixscape": [
        ("perturbation_signature.allclose_excess", "<=", 1.0, ""),
        ("perturbation_signature.pearson_correlation", ">=", 0.999, ""),
        ("mixscape.global_class.exact_agreement", ">=", 0.95, ""),
        ("mixscape.p_ko.pearson_correlation", ">=", 0.95, ""),
        ("mixscape.lda.lda_abs_correlation", ">=", 0.95, ""),
    ],
    "rank_genes_groups": [
        ("*.score.pearson_correlation", ">=", 0.98, ""),
        ("*.markers.set_jaccard", ">=", 0.85, ""),
    ],
    "scanpy_core_graphs_embeddings": [
        ("neighbors.distance.graph_jaccard", ">=", 0.9, ""),
        ("neighbors.connectivity.graph_jaccard", ">=", 0.85, ""),
        ("umap.cross_embedding_knn_overlap", ">=", 0.65, ""),
        ("leiden.adjusted_rand_index", ">=", 0.9, ""),
        ("leiden.normalized_mutual_information", ">=", 0.9, ""),
    ],
    "scanpy_core_harmony": [
        ("harmony.minimum_component_abs_correlation", ">=", 0.95, ALLCLOSE_BASIS),
        ("harmony.mean_component_abs_correlation", ">=", 0.98, ""),
        ("harmony.standard_deviation_max_abs_error", "<=", 0.1, ""),
    ],
    "scanpy_core_hvg_pca": [
        ("highly_variable_genes.seurat.selection.set_jaccard", ">=", 0.99, ""),
        ("highly_variable_genes.seurat.means.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.seurat.dispersions.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.seurat.dispersions_norm.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.cell_ranger.selection.set_jaccard", ">=", 0.99, ""),
        ("highly_variable_genes.cell_ranger.means.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.cell_ranger.dispersions.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.cell_ranger.dispersions_norm.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.seurat_v3.selection.set_jaccard", ">=", 0.99, ""),
        ("highly_variable_genes.seurat_v3.means.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.seurat_v3.variances.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.seurat_v3.variances_norm.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.pearson_residuals.selection.set_jaccard", ">=", 0.99, ""),
        ("highly_variable_genes.pearson_residuals.means.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.pearson_residuals.variances.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("highly_variable_genes.pearson_residuals.residual_variances.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
        ("pca.scores.minimum_component_abs_correlation", ">=", 0.95, ""),
        ("pca.loadings.minimum_component_abs_correlation", ">=", 0.95, ""),
        ("pca.explained_variance_ratio.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
    ],
    "scanpy_core_preprocessing": [
        ("filter_cells.index_agreement", ">=", 1.0, ""),
        ("filter_genes.index_agreement", ">=", 1.0, ""),
        ("calculate_qc_metrics.*.allclose_excess", "<=", 1.0, ALLCLOSE_BASIS),
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


# Why a criterion fails, where that has been investigated. Kept here rather than in the
# comparison scripts because it is a statement about the criterion, not a measurement, and
# because `evaluate.py` shows it only for metrics that actually fail.
DIAGNOSES: dict[tuple[str, str], str] = {
    ("biological_pipeline_pbmc3k", "umap.cross_embedding_knn_overlap"): (
        "UMAP is stochastic, and this threshold asks for more agreement than the CPU reference shows "
        "against itself. Compare `umap.cpu_reseeded_knn_overlap`, the same measurement with only the "
        "seed changed, and `umap.cross_embedding_overlap_vs_cpu_baseline`."
    ),
    ("scanpy_core_graphs_embeddings", "umap.cross_embedding_knn_overlap"): (
        "UMAP is stochastic, and this threshold asks for more agreement than the CPU reference shows "
        "against itself. On this dataset the GPU embedding is closer to the CPU one than a reseeded "
        "CPU run is, and the criterion still fails."
    ),
    ("calculate_niche", "neighborhood.adjusted_rand_index"): (
        "The neighborhood profile holds only 755 distinct rows across 4668 cells, so 89.67% of cells "
        "have an exact distance tie at the k-th neighbour. Against exact float64 ground truth the GPU "
        "kNN is exact and `sc.pp.neighbors` is not, on 743 rows: the divergence is the CPU reference."
    ),
    ("calculate_niche", "utag.adjusted_rand_index"): (
        "Squidpy calls `sc.tl.leiden` without a flavor, so CPU and GPU use different Leiden backends. "
        "On identical input Scanpy's own leidenalg and igraph backends agree only at ARI 0.5041, below "
        "this threshold, so the criterion measures backend choice rather than correctness."
    ),
    ("calculate_niche", "*.cluster_count_difference"): (
        "`resolution` does not carry the same meaning across Leiden implementations: on identical input "
        "cuGraph found 34 clusters where leidenalg found 41, and Scanpy's own two backends already "
        "differ by 1."
    ),
    ("*", "*.allclose_excess"): (
        "numpy.allclose has two terms, atol + rtol * |b|, and which one binds depends on the magnitude "
        "of the element that fails. Where a quantity passes through zero the absolute floor binds, and "
        "at float32 precision no implementation can satisfy it there. See NUMERICAL_VALIDATION.md for "
        "the per-comparison split between that case and a genuine relative disagreement."
    ),
}
DIAGNOSES[("calculate_niche", "neighborhood.normalized_mutual_information")] = DIAGNOSES[
    ("calculate_niche", "neighborhood.adjusted_rand_index")
]
DIAGNOSES[("calculate_niche", "utag.normalized_mutual_information")] = DIAGNOSES[
    ("calculate_niche", "utag.adjusted_rand_index")
]


def diagnosis_for(method: str, metric: str) -> str:
    """The recorded explanation for a failing criterion, or "" if none was written."""
    for key in ((method, metric), ("*", metric)):
        if key in DIAGNOSES:
            return DIAGNOSES[key]
    for (rule_method, pattern), text in DIAGNOSES.items():
        if rule_method in (method, "*") and fnmatch.fnmatchcase(metric, pattern):
            return text
    return ""
