# Decoupler equivalence comparisons

`decoupler.py` runs MLM, ULM, AUCell, WAGGR, and Z-score against the matching
`decoupler.mt` methods on the same seeded toy data and network. It compares output
schemas, scores, and adjusted p-values where the method produces them.

WAGGR uses `times=0` here to isolate deterministic weighted aggregation from
permutation sampling. Run the script in the project GPU environment; its JSON report
is written to `results/` or `$EQUIVALENCE_OUTPUT_DIR`.
