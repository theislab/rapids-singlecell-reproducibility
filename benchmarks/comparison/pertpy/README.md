# Pertpy equivalence comparisons

The scripts compare every documented rapids-singlecell Pertpy-equivalent API group
against `pertpy`: all supported `Distance` metrics, three `GuideAssignment` methods,
Mixscape perturbation signatures/classification/LDA, and Mixscale scores.

Mixscale requires a Pertpy release that exposes `pt.tl.Mixscale`; the script fails with
an actionable message instead of silently treating a missing CPU reference as success.
Run these scripts in the project GPU environment. They store both implementations' outputs
for `../evaluate.py`; records go to `results/` or `$EQUIVALENCE_OUTPUT_DIR`.
