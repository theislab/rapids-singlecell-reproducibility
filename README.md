# rapids_singlecell-reproducibility

This repository contains all scripts to reproduce the associated publication.

## Use cases

TODO

## Reproducibility

Method-level CPU/GPU comparisons live in [`benchmarks/comparison`](benchmarks/comparison).
The original Scanpy comparisons are complemented by structured comparisons for the
previously uncovered Scanpy, Squidpy, Decoupler, and Pertpy APIs:

- [`scanpy_extended`](benchmarks/comparison/scanpy_extended/README.md)
- [`squidpy`](benchmarks/comparison/squidpy/README.md)
- [`decoupler`](benchmarks/comparison/decoupler/README.md)
- [`pertpy`](benchmarks/comparison/pertpy/README.md)

Each new script writes a JSON record containing the method, reference-package version,
dataset, metric, tolerance, observed value, tier, and pass/fail status. After running
the comparisons individually, aggregate the records with:

```bash
python benchmarks/comparison/collect_results.py
```

To run the complete structured suite in isolated processes and aggregate it
automatically:

```bash
python benchmarks/comparison/run_structured.py
```

For the Theislab/HMGU GPU cluster, use the documented
[SLURM + uv + localscratch workflow](cluster/README.md).

See the complete [CPU/GPU coverage inventory](benchmarks/comparison/COVERAGE.md) for
the mapping from public methods to evidence scripts.

## Run time

TODO
