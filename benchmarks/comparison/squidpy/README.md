# Squidpy equivalence comparisons

These scripts compare the accelerated `rapids-singlecell` implementations directly
against the public `squidpy.gr` implementations. Each script stores both implementations' outputs and
writes a record to `results/` (or `$EQUIVALENCE_OUTPUT_DIR`); `../evaluate.py` compares them
and reports a missed criterion.

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
OmniPath download. The deterministic methods are scored with `numpy.allclose` at its default parameters and a
correlation on the same arrays.
For `ligrec`, deterministic means are checked tightly while independently sampled
permutation p-values use distribution-level tolerances.

`calculate_niche` compares all three flavors implemented by rapids-singlecell:
neighborhood, UTAG, and CellCharter. Since cluster identifiers are arbitrary, it uses
label-invariant ARI/NMI plus cluster-count agreement.
