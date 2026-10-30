# CPU/GPU equivalence report

Scored 2026-10-30.

## Outcome

- Overall: **PASS**
- Gating metrics passing: **161/161**
- Scripts completing successfully: **17/17**

## Software versions

| Package | Version(s) |
| --- | --- |
| decoupler | 2.2.0 |
| harmonypy | 2.0.2 |
| pertpy | 1.3.0 |
| rapids-singlecell | 0.17.0 |
| scanpy | 1.13.0a2 |
| squidpy | 1.8.4 |

## Hardware

| | |
| --- | --- |
| GPU | NVIDIA A100 80GB PCIe |
| Compute capability | 8.0 |
| CUDA driver / runtime | 13.3 / 13.2 |
| Device memory | 79.3 GiB |
| Measured | 2026-10-30T09:01:15.083859+00:00 |

## Method groups

| Method group | Reference | Dataset | Scale | Result | Metrics |
| --- | --- | --- | --- | --- | ---: |
| `bbknn_scrublet` | scanpy | pbmc68k_reduced + pbmc3k |  | PASS | 5/5 |
| `biological_pipeline_pbmc3k` | scanpy | pbmc3k_processed raw log-expression with published cell-type labels | 2,638 x 13,714 | PASS | 10/10 |
| `co_occurrence` | squidpy | squidpy.datasets.imc | 4,668 x 34 | PASS | 2/2 |
| `decoupler_methods` | decoupler | decoupler.ds.toy | 80 x 40 | PASS | 9/9 |
| `distance` | pertpy | seeded grouped Gaussian data | 120 x 12 | PASS | 18/18 |
| `embeddings_extended` | scanpy | scanpy.datasets.pbmc68k_reduced | 700 x 765 | PASS | 6/6 |
| `guide_assignment` | pertpy | seeded Poisson guide-count mixture | 200 x 8 | PASS | 3/3 |
| `ingest_cell_cycle` | scanpy | scanpy.datasets.pbmc68k_reduced | 700 x 765 | PASS | 6/6 |
| `ligrec` | squidpy | scanpy.datasets.paul15 | 742 x 3,451 | PASS | 4/4 |
| `mixscale` | pertpy | seeded synthetic perturbation screen | 30 x 20 | PASS | 1/1 |
| `mixscape` | pertpy | seeded synthetic perturbation screen | 30 x 20 | PASS | 4/4 |
| `rank_genes_groups` | scanpy | scanpy.datasets.pbmc68k_reduced | 700 x 765 | PASS | 40/40 |
| `scanpy_core_graphs_embeddings` | scanpy | pbmc68k_reduced | 700 x 765 | PASS | 5/5 |
| `scanpy_core_harmony` | harmonypy | Harmonypy PBMC 3,500-cell donor benchmark | 3,500 x 30 | PASS | 4/4 |
| `scanpy_core_hvg_pca` | scanpy | pbmc3k | 2,700 x 13,714 | PASS | 21/21 |
| `scanpy_core_preprocessing` | scanpy | pbmc3k | 2,700 x 32,738 | PASS | 21/21 |
| `spatial_autocorr` | squidpy | squidpy.datasets.imc | 4,668 x 34 | PASS | 2/2 |

## Criteria

| Method group | Metric | Observed | Criterion | Result |
| --- | --- | ---: | --- | --- |
| `bbknn_scrublet` | `bbknn.graph_jaccard` | 1 | >= 0.85 | PASS |
| `bbknn_scrublet` | `scrublet.call.exact_agreement` | 0.9962963 | >= 0.95 | PASS |
| `bbknn_scrublet` | `scrublet.score.pearson_correlation` | 0.97489925 | >= 0.95 | PASS |
| `bbknn_scrublet` | `scrublet_simulate_doublets.library_size.relative_error_of_mean` | 0 | <= 0.02 | PASS |
| `bbknn_scrublet` | `scrublet_simulate_doublets.library_size.relative_error_of_std` | 0 | <= 0.05 | PASS |
| `biological_pipeline_pbmc3k` | `umap.cross_embedding_overlap_vs_cpu_baseline` | 0.020217336 | >= -0.05 | PASS |
| `biological_pipeline_pbmc3k` | `annotation.exact_agreement` | 1 | >= 0.9 | PASS |
| `biological_pipeline_pbmc3k` | `annotation.accuracy_difference` | 0 | <= 0.05 | PASS |
| `biological_pipeline_pbmc3k` | `clustering.adjusted_rand_index` | 0.97825978 | >= 0.8 | PASS |
| `biological_pipeline_pbmc3k` | `clustering.normalized_mutual_information` | 0.96283526 | >= 0.8 | PASS |
| `biological_pipeline_pbmc3k` | `clustering.cell_type_nmi_difference` | 0.0051196805 | <= 0.05 | PASS |
| `biological_pipeline_pbmc3k` | `highly_variable_genes.selection.set_jaccard` | 1 | >= 0.98 | PASS |
| `biological_pipeline_pbmc3k` | `pca.minimum_component_abs_correlation` | 0.99999927 | >= 0.95 | PASS |
| `biological_pipeline_pbmc3k` | `markers.mean_set_jaccard` | 1 | >= 0.85 | PASS |
| `biological_pipeline_pbmc3k` | `markers.minimum_set_jaccard` | 1 | >= 0.7 | PASS |
| `co_occurrence` | `occurrence.mean_abs_error` | 9.8388647e-08 | <= 1e-06 | PASS |
| `co_occurrence` | `occurrence.pearson_correlation` | 1 | >= 0.99999 | PASS |
| `decoupler_methods` | `aucell.score.pearson_correlation` | 1 | >= 0.999 | PASS |
| `decoupler_methods` | `mlm.adjusted_pvalue.pearson_correlation` | 1 | >= 0.999 | PASS |
| `decoupler_methods` | `mlm.score.pearson_correlation` | 1 | >= 0.999 | PASS |
| `decoupler_methods` | `ulm.adjusted_pvalue.pearson_correlation` | 1 | >= 0.999 | PASS |
| `decoupler_methods` | `ulm.score.pearson_correlation` | 1 | >= 0.999 | PASS |
| `decoupler_methods` | `waggr.adjusted_pvalue.pearson_correlation` | 1 | >= 0.999 | PASS |
| `decoupler_methods` | `waggr.score.pearson_correlation` | 1 | >= 0.999 | PASS |
| `decoupler_methods` | `zscore.adjusted_pvalue.pearson_correlation` | 1 | >= 0.999 | PASS |
| `decoupler_methods` | `zscore.score.pearson_correlation` | 1 | >= 0.999 | PASS |
| `distance` | `cosine_distance.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `cosine_distance.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `edistance.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `edistance.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `euclidean.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `euclidean.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `mean_absolute_error.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `mean_absolute_error.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `mse.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `mse.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `pearson_distance.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `pearson_distance.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `r2_distance.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `r2_distance.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `root_mean_squared_error.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `root_mean_squared_error.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `wasserstein.onesided.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `distance` | `wasserstein.pairwise.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `embeddings_extended` | `diffmap.minimum_component_abs_correlation` | 1 | >= 0.95 | PASS |
| `embeddings_extended` | `draw_graph.cross_embedding_knn_overlap` | 0.61371429 | >= 0.6 | PASS |
| `embeddings_extended` | `embedding_density.pearson_correlation` | 1 | >= 0.999 | PASS |
| `embeddings_extended` | `tsne.cross_embedding_knn_overlap` | 0.48504762 | >= 0.45 | PASS |
| `embeddings_extended` | `tsne.cpu.trustworthiness` | 0.94223113 | >= 0.9 | PASS |
| `embeddings_extended` | `tsne.gpu.trustworthiness` | 0.92014827 | >= 0.9 | PASS |
| `guide_assignment` | `assign_by_threshold.exact_agreement` | 1 | >= 1.0 | PASS |
| `guide_assignment` | `assign_mixture_model.exact_agreement` | 1 | >= 0.9 | PASS |
| `guide_assignment` | `assign_to_max_guide.exact_agreement` | 1 | >= 1.0 | PASS |
| `ingest_cell_cycle` | `ingest.label.exact_agreement` | 0.99142857 | >= 0.98 | PASS |
| `ingest_cell_cycle` | `ingest.pca.pearson_correlation` | 1 | >= 0.999 | PASS |
| `ingest_cell_cycle` | `ingest.umap.cross_embedding_knn_overlap` | 0.50095238 | >= 0.45 | PASS |
| `ingest_cell_cycle` | `score_genes_cell_cycle.G2M.pearson_correlation` | 1 | >= 0.999 | PASS |
| `ingest_cell_cycle` | `score_genes_cell_cycle.S.pearson_correlation` | 1 | >= 0.999 | PASS |
| `ingest_cell_cycle` | `score_genes_cell_cycle.phase.exact_agreement` | 1 | >= 0.99 | PASS |
| `ligrec` | `means.float32_allclose_excess` | 0.027717381 | <= 1.0 | PASS |
| `ligrec` | `pvalues.mean_abs_error` | 0.015612246 | <= 0.05 | PASS |
| `ligrec` | `pvalues.pearson_correlation` | 0.99810767 | >= 0.9 | PASS |
| `ligrec` | `pvalues.nan_mask_agreement` | 1 | >= 1.0 | PASS |
| `mixscale` | `mixscale.score.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `mixscape` | `mixscape.global_class.exact_agreement` | 1 | >= 0.95 | PASS |
| `mixscape` | `mixscape.lda.lda_abs_correlation` | 1 | >= 0.95 | PASS |
| `mixscape` | `mixscape.p_ko.pearson_correlation` | 1 | >= 0.95 | PASS |
| `mixscape` | `perturbation_signature.pearson_correlation` | 1 | >= 0.999 | PASS |
| `rank_genes_groups` | `t-test.CD14+ Monocyte.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD14+ Monocyte.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.CD19+ B.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD19+ B.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.CD34+.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD34+.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.CD4+_CD25 T Reg.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD4+_CD25 T Reg.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.CD4+_CD45RA+_CD25- Naive T.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD4+_CD45RA+_CD25- Naive T.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.CD4+_CD45RO+ Memory.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD4+_CD45RO+ Memory.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.CD56+ NK.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD56+ NK.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.CD8+ Cytotoxic T.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD8+ Cytotoxic T.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.CD8+_CD45RA+ Naive Cytotoxic.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.CD8+_CD45RA+ Naive Cytotoxic.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `t-test.Dendritic.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `t-test.Dendritic.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD14+ Monocyte.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD14+ Monocyte.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD19+ B.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD19+ B.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD34+.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD34+.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD4+_CD25 T Reg.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD4+_CD25 T Reg.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD4+_CD45RA+_CD25- Naive T.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD4+_CD45RA+_CD25- Naive T.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD4+_CD45RO+ Memory.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD4+_CD45RO+ Memory.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD56+ NK.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD56+ NK.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD8+ Cytotoxic T.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD8+ Cytotoxic T.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.CD8+_CD45RA+ Naive Cytotoxic.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.CD8+_CD45RA+ Naive Cytotoxic.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `rank_genes_groups` | `wilcoxon.Dendritic.markers.set_jaccard` | 1 | >= 0.85 | PASS |
| `rank_genes_groups` | `wilcoxon.Dendritic.score.pearson_correlation` | 1 | >= 0.98 | PASS |
| `scanpy_core_graphs_embeddings` | `umap.cross_embedding_overlap_vs_cpu_baseline` | 0.032190476 | >= -0.05 | PASS |
| `scanpy_core_graphs_embeddings` | `leiden.adjusted_rand_index` | 0.99130605 | >= 0.9 | PASS |
| `scanpy_core_graphs_embeddings` | `leiden.normalized_mutual_information` | 0.98796652 | >= 0.9 | PASS |
| `scanpy_core_graphs_embeddings` | `neighbors.connectivity.graph_jaccard` | 1 | >= 0.85 | PASS |
| `scanpy_core_graphs_embeddings` | `neighbors.distance.graph_jaccard` | 1 | >= 0.9 | PASS |
| `scanpy_core_harmony` | `harmony.standard_deviation_max_abs_error` | 0.024928254 | <= 0.1 | PASS |
| `scanpy_core_harmony` | `harmony.relative_l2_max` | 0.070319773 | <= 0.1 | PASS |
| `scanpy_core_harmony` | `harmony.minimum_component_abs_correlation` | 0.99761202 | > 0.95 | PASS |
| `scanpy_core_harmony` | `harmony.mean_component_abs_correlation` | 0.99886787 | >= 0.98 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.cell_ranger.dispersions.float32_allclose_excess` | 0.0043552499 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.cell_ranger.dispersions_norm.float32_allclose_excess` | 0.3685738 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.cell_ranger.means.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.cell_ranger.selection.set_jaccard` | 1 | >= 0.99 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.pearson_residuals.means.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.pearson_residuals.residual_variances.float32_allclose_excess` | 0.33199474 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.pearson_residuals.selection.set_jaccard` | 1 | >= 0.99 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.pearson_residuals.variances.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat.dispersions.float32_allclose_excess` | 0.0093441301 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat.dispersions_norm.float32_allclose_excess` | 0.098022152 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat.means.float32_allclose_excess` | 0.002614238 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat.selection.set_jaccard` | 1 | >= 0.99 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat_v3.means.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat_v3.selection.set_jaccard` | 1 | >= 0.99 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat_v3.variances.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat_v3.variances_norm.float32_allclose_excess` | 1.515351e-08 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `pca.explained_variance_ratio.float32_allclose_excess` | 4.384644e-11 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `pca.loadings.abs_float32_allclose_excess` | 2.5812364e-07 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `pca.loadings.minimum_component_abs_correlation` | 1 | >= 0.95 | PASS |
| `scanpy_core_hvg_pca` | `pca.scores.abs_float32_allclose_excess` | 0.0071307612 | <= 1.0 | PASS |
| `scanpy_core_hvg_pca` | `pca.scores.minimum_component_abs_correlation` | 1 | >= 0.95 | PASS |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.n_genes_by_counts.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.pct_counts_hb.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.pct_counts_mt.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.pct_counts_ribo.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.total_counts.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.total_counts_hb.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.total_counts_mt.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.total_counts_ribo.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `filter_cells.index_agreement` | 1 | >= 1.0 | PASS |
| `scanpy_core_preprocessing` | `filter_genes.index_agreement` | 1 | >= 1.0 | PASS |
| `scanpy_core_preprocessing` | `log1p.pearson_correlation` | 1 | >= 0.999999 | PASS |
| `scanpy_core_preprocessing` | `log1p.float32_allclose_excess` | 0.013653004 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `normalize_pearson_residuals.pearson_correlation` | 1 | >= 0.99999 | PASS |
| `scanpy_core_preprocessing` | `normalize_pearson_residuals.float32_allclose_excess` | 0.52745774 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `normalize_total.pearson_correlation` | 1 | >= 0.999999 | PASS |
| `scanpy_core_preprocessing` | `normalize_total.float32_allclose_excess` | 0.014913623 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `regress_out.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `scanpy_core_preprocessing` | `scale.pearson_correlation` | 1 | >= 0.99999 | PASS |
| `scanpy_core_preprocessing` | `scale.float32_allclose_excess` | 0.032224136 | <= 1.0 | PASS |
| `scanpy_core_preprocessing` | `score_genes.pearson_correlation` | 1 | >= 0.9999 | PASS |
| `scanpy_core_preprocessing` | `sqrt.float32_allclose_excess` | 0 | <= 1.0 | PASS |
| `spatial_autocorr` | `geary.C.pearson_correlation` | 1 | >= 0.999999 | PASS |
| `spatial_autocorr` | `moran.I.pearson_correlation` | 1 | >= 0.999999 | PASS |

## Recorded measurements (not gating)

| Method group | Measurement | Observed |
| --- | --- | ---: |
| `biological_pipeline_pbmc3k` | `umap.cpu_reseeded_knn_overlap` | 0.40328532 |
| `biological_pipeline_pbmc3k` | `markers.B cells.set_jaccard` | 1 |
| `biological_pipeline_pbmc3k` | `markers.CD14+ Monocytes.set_jaccard` | 1 |
| `biological_pipeline_pbmc3k` | `markers.CD4 T cells.set_jaccard` | 1 |
| `biological_pipeline_pbmc3k` | `markers.CD8 T cells.set_jaccard` | 1 |
| `biological_pipeline_pbmc3k` | `markers.Dendritic cells.set_jaccard` | 1 |
| `biological_pipeline_pbmc3k` | `markers.FCGR3A+ Monocytes.set_jaccard` | 1 |
| `biological_pipeline_pbmc3k` | `markers.Megakaryocytes.set_jaccard` | 1 |
| `biological_pipeline_pbmc3k` | `markers.NK cells.set_jaccard` | 1 |
| `biological_pipeline_pbmc3k` | `umap.cross_embedding_knn_overlap` | 0.42350265 |
| `co_occurrence` | `occurrence.allclose_excess` | 0.047047078 |
| `decoupler_methods` | `aucell.score.allclose_excess` | 0.0051959217 |
| `decoupler_methods` | `mlm.adjusted_pvalue.allclose_excess` | 1.820935 |
| `decoupler_methods` | `mlm.score.allclose_excess` | 10.356613 |
| `decoupler_methods` | `ulm.adjusted_pvalue.allclose_excess` | 21.178564 |
| `decoupler_methods` | `ulm.score.allclose_excess` | 0.62557934 |
| `decoupler_methods` | `waggr.adjusted_pvalue.allclose_excess` | 0 |
| `decoupler_methods` | `waggr.score.allclose_excess` | 0.71895189 |
| `decoupler_methods` | `zscore.adjusted_pvalue.allclose_excess` | 0.14219173 |
| `decoupler_methods` | `zscore.score.allclose_excess` | 0.053337029 |
| `distance` | `cosine_distance.onesided.allclose_excess` | 2.1415124e-11 |
| `distance` | `cosine_distance.pairwise.allclose_excess` | 5.5241235e-10 |
| `distance` | `edistance.onesided.allclose_excess` | 2.5593545e-09 |
| `distance` | `edistance.pairwise.allclose_excess` | 1.020402e-09 |
| `distance` | `euclidean.onesided.allclose_excess` | 1.7046396e-11 |
| `distance` | `euclidean.pairwise.allclose_excess` | 4.8616687e-11 |
| `distance` | `mean_absolute_error.onesided.allclose_excess` | 2.1716695e-11 |
| `distance` | `mean_absolute_error.pairwise.allclose_excess` | 2.2845928e-11 |
| `distance` | `mse.onesided.allclose_excess` | 2.3046615e-11 |
| `distance` | `mse.pairwise.allclose_excess` | 1.1945331e-10 |
| `distance` | `pearson_distance.onesided.allclose_excess` | 6.9996527e-11 |
| `distance` | `pearson_distance.pairwise.allclose_excess` | 1.2249392e-10 |
| `distance` | `r2_distance.onesided.allclose_excess` | 1.4076397e-10 |
| `distance` | `r2_distance.pairwise.allclose_excess` | 2.3490795e-10 |
| `distance` | `root_mean_squared_error.onesided.allclose_excess` | 1.7046396e-11 |
| `distance` | `root_mean_squared_error.pairwise.allclose_excess` | 1.5800423e-10 |
| `distance` | `wasserstein.onesided.allclose_excess` | 0.42112604 |
| `distance` | `wasserstein.pairwise.allclose_excess` | 0.9864693 |
| `ligrec` | `means.allclose_excess` | 0.036598376 |
| `mixscale` | `mixscale.score.allclose_excess` | 0.0051948259 |
| `mixscape` | `perturbation_signature.allclose_excess` | 2.6862685 |
| `scanpy_core_graphs_embeddings` | `umap.cpu_reseeded_knn_overlap` | 0.55466667 |
| `scanpy_core_graphs_embeddings` | `neighbors.connectivity.graph_exact_agreement` | 1 |
| `scanpy_core_graphs_embeddings` | `neighbors.distance.graph_exact_agreement` | 1 |
| `scanpy_core_graphs_embeddings` | `umap.cross_embedding_knn_overlap` | 0.58685714 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.cell_ranger.dispersions.allclose_excess` | 0.0048088648 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.cell_ranger.dispersions_norm.allclose_excess` | 13.067911 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.cell_ranger.means.allclose_excess` | 0 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.pearson_residuals.means.allclose_excess` | 0 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.pearson_residuals.residual_variances.allclose_excess` | 0.35402906 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.pearson_residuals.variances.allclose_excess` | 0 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat.dispersions.allclose_excess` | 0.037546118 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat.dispersions_norm.allclose_excess` | 3.8232139 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat.means.allclose_excess` | 0.0067183445 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat_v3.means.allclose_excess` | 0 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat_v3.variances.allclose_excess` | 0 |
| `scanpy_core_hvg_pca` | `highly_variable_genes.seurat_v3.variances_norm.allclose_excess` | 1.606809e-08 |
| `scanpy_core_hvg_pca` | `pca.explained_variance_ratio.allclose_excess` | 4.2991796e-10 |
| `scanpy_core_hvg_pca` | `pca.loadings.abs_allclose_excess` | 1.0342365e-05 |
| `scanpy_core_hvg_pca` | `pca.scores.abs_allclose_excess` | 0.0094870597 |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.n_genes_by_counts.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.pct_counts_hb.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.pct_counts_mt.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.pct_counts_ribo.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.total_counts.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.total_counts_hb.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.total_counts_mt.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `calculate_qc_metrics.total_counts_ribo.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `log1p.allclose_excess` | 0.014473541 |
| `scanpy_core_preprocessing` | `normalize_pearson_residuals.allclose_excess` | 30.626773 |
| `scanpy_core_preprocessing` | `normalize_total.allclose_excess` | 0.014920842 |
| `scanpy_core_preprocessing` | `regress_out.allclose_excess` | 12.018964 |
| `scanpy_core_preprocessing` | `scale.allclose_excess` | 1.7870379 |
| `scanpy_core_preprocessing` | `score_genes.allclose_excess` | 0 |
| `scanpy_core_preprocessing` | `sqrt.allclose_excess` | 0 |
| `spatial_autocorr` | `geary.C.allclose_excess` | 0.0013172383 |
| `spatial_autocorr` | `moran.I.allclose_excess` | 0.0011057476 |
