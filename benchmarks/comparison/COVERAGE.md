# CPU/GPU equivalence coverage

This inventory maps the rapids-singlecell API surface to manuscript-facing CPU/GPU
comparisons. “Existing” denotes the original one-method scripts; “structured” denotes
new scripts that emit JSON records with explicit metrics and thresholds.

| API group | Covered methods | Evidence |
|---|---|---|
| Scanpy preprocessing | `filter_cells`, `filter_genes`, `calculate_qc_metrics`, `normalize_total`, `log1p`, `normalize_pearson_residuals`, all four HVG flavors, `regress_out`, `scale`, `pca`, `harmony_integrate`, `neighbors`, `score_genes` | `scanpy_core` structured scripts; original one-method assertions retained |
| Additional Scanpy preprocessing | `sqrt`, `bbknn`, `scrublet`, `scrublet_simulate_doublets` | `scanpy_extended` structured scripts |
| Scanpy tools | `leiden`, `umap` | Existing scripts, strengthened with label/embedding agreement metrics |
| Additional Scanpy tools | `louvain`, `kmeans`, `tsne`, `diffmap`, `draw_graph`, `embedding_density`, `ingest`, `rank_genes_groups` (t-test, Wilcoxon, logistic regression), `score_genes_cell_cycle` | `scanpy_extended` structured scripts |
| Squidpy | `spatial_autocorr` (Moran and Geary), `co_occurrence`, `ligrec`, `calculate_niche` (neighborhood, UTAG, CellCharter) | `squidpy` structured scripts |
| Decoupler | `mlm`, `ulm`, `aucell`, `waggr`, `zscore` | `decoupler` structured script |
| Pertpy Distance | All nine supported metrics through `pairwise` | `pertpy/distance.py` |
| Pertpy GuideAssignment | Threshold, maximum-guide, and mixture-model assignment | `pertpy/guide_assignment.py` |
| Pertpy Mixscape | Perturbation signature, classification, and LDA | `pertpy/mixscape.py` |
| Pertpy Mixscale | Perturbation signature (shared with Mixscape) and continuous Mixscale score | `pertpy/mixscale.py` |
| End-to-end biological interpretation | Independent HVG, scaling, PCA, neighbors, UMAP, Leiden, marker ranking, and held-out annotation workflows | `biological_pipeline/pbmc3k.py` |

`kmeans` is included because it is exported by `rapids_singlecell.tl`, although it is
not currently listed in `docs/api/scanpy_gpu.md`. Squidpy's `spatialleiden` flavor is
not listed because rapids-singlecell does not implement it.

The structured reports can be aggregated into `equivalence.json` with
`python benchmarks/comparison/collect_results.py` after running them in the GPU
environment.

`python benchmarks/comparison/run_structured.py` executes the complete suite in isolated
processes, retains all failures, and produces a Markdown report, flat metric table, pass-rate
figure, biological UMAP comparison, and marker-overlap table. The expensive
`calculate_niche` comparison runs last.
