# End-to-end biological equivalence

`pbmc3k.py` starts from the published PBMC3k log-expression matrix and cell-type labels,
then runs independent CPU and GPU workflows through HVG selection, scaling, PCA, neighbor
graph construction, UMAP, Leiden clustering, and marker ranking.

The comparison reports:

- HVG overlap and PCA component agreement;
- UMAP trustworthiness and cross-embedding neighborhood preservation;
- CPU/GPU clustering ARI and NMI, plus agreement with the published cell types;
- mean and worst-case top-50 marker overlap per cell type; and
- held-out kNN annotation accuracy, CPU/GPU prediction agreement, and accuracy difference.

It also writes a side-by-side UMAP figure and per-cell-type marker-overlap CSV to
`$EQUIVALENCE_ARTIFACT_DIR`. This is the manuscript-facing biological interpretation test;
the smaller method tests remain useful for diagnosing where a pipeline divergence begins.
