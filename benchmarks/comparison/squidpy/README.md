# Squidpy equivalence comparisons

These scripts compare the accelerated `rapids-singlecell` implementations directly
against the public `squidpy.gr` implementations. Each script writes a structured JSON
record to `results/` (or `$EQUIVALENCE_OUTPUT_DIR`) and exits non-zero when an explicit
equivalence criterion is missed.

Run them in a Linux environment with an NVIDIA GPU and `rapids-singlecell`, `squidpy`,
`scanpy`, and their dependencies installed:

```bash
python benchmarks/comparison/squidpy/spatial_autocorr.py
python benchmarks/comparison/squidpy/co_occurrence.py
python benchmarks/comparison/squidpy/ligrec.py
python benchmarks/comparison/squidpy/calculate_niche.py
```

`spatial_autocorr` and `co_occurrence` use Squidpy's labeled IMC example dataset.
`ligrec` uses the Paul15 dataset and a fixed local interaction set, avoiding a mutable
OmniPath download. The deterministic methods use direct numerical error and correlation.
For `ligrec`, deterministic means are checked tightly while independently sampled
permutation p-values use distribution-level tolerances.

`calculate_niche` compares all three flavors implemented by rapids-singlecell:
neighborhood, UTAG, and CellCharter. Since cluster identifiers are arbitrary, it uses
label-invariant ARI/NMI plus cluster-count agreement.
