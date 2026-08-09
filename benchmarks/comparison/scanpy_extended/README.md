# Additional Scanpy equivalence comparisons

These scripts cover the previously missing Scanpy-equivalent methods: BBKNN,
Scrublet (including doublet simulation), Louvain, k-means, t-SNE, DiffMap,
force-directed graph drawing, embedding density, ingest, cell-cycle scoring, and
marker ranking with t-test, Wilcoxon, and logistic regression.

Run each script in the project GPU environment. Each stores the CPU and GPU outputs to
`$EQUIVALENCE_ARRAY_DIR` and its record to `results/` or `$EQUIVALENCE_OUTPUT_DIR`;
`../evaluate.py` does the comparing.
