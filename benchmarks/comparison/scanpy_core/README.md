# Core Scanpy equivalence comparisons

These scripts replace the original assertion-only checks with reviewer-facing structured
records for the primary methods discussed in the manuscript:

- `preprocessing.py`: filtering, QC, normalization, log transformation, Pearson residuals,
  scaling, regression, gene scoring, and square root transformation;
- `hvg_pca.py`: all four supported HVG flavors and sign-invariant PCA agreement;
- `graphs_embeddings.py`: neighbor graph overlap, UMAP fidelity, and Leiden ARI/NMI; and
- `harmony.py`: a live Harmonypy 0.2.0 CPU run on Harmonypy's 3,500-cell donor
  benchmark against the GPU implementation's original-Harmony flavor. It calls
  Harmonypy directly because Scanpy 1.12's wrapper assumes the output orientation of
  Harmonypy 0.0.x.

Every script writes observed metrics and explicit thresholds to `$EQUIVALENCE_OUTPUT_DIR`
before raising for a failed threshold. Run them through `../run_structured.py` so failures,
logs, software versions, and manuscript-ready tables are collected consistently.
