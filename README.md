# rapids_singlecell-reproducibility

This repository contains all scripts to reproduce the associated publication.

## Use cases

- [`use_cases/100M-working.ipynb`](use_cases/100M-working.ipynb) contains the Dask-based
  Tahoe-100M workflow.
- [`benchmarks/speed/notebooks`](benchmarks/speed/notebooks) contains matching CPU and GPU
  notebooks for the million-cell mouse-brain workflow.

## Reproducibility

Method-level CPU/GPU comparisons live in [`benchmarks/comparison`](benchmarks/comparison).
The original Scanpy comparisons are complemented by structured comparisons for the
previously uncovered Scanpy, Squidpy, Decoupler, and Pertpy APIs:

- [`scanpy_extended`](benchmarks/comparison/scanpy_extended/README.md)
- [`squidpy`](benchmarks/comparison/squidpy/README.md)
- [`decoupler`](benchmarks/comparison/decoupler/README.md)
- [`pertpy`](benchmarks/comparison/pertpy/README.md)

Each new script writes a JSON record containing the method, reference-package version,
dataset, metric, tolerance, observed value, tier, and pass/fail status. Metrics are split
into criteria that gate the suite and measurements that are recorded as evidence. Thresholds are
never widened to make a run green; failing criteria are diagnosed and left in place. See
[`THRESHOLDS.md`](benchmarks/comparison/THRESHOLDS.md) for the record schema and for the
diagnosis behind each current failure. After running the comparisons individually, aggregate the
records with:

```bash
python benchmarks/comparison/collect_results.py
```

To run the complete structured suite in isolated processes and aggregate it
automatically:

```bash
python benchmarks/comparison/run_structured.py
```

Pass a subset of the script paths to rerun only those comparisons, which marks the aggregate
output as partial:

```bash
python benchmarks/comparison/run_structured.py scanpy_core/preprocessing.py
```

The complete suite includes deterministic numerical comparisons, stochastic graph and
embedding comparisons, the accelerated Squidpy/Decoupler/Pertpy APIs, and an end-to-end
PBMC3k biological workflow. It always preserves per-script logs and produces:

- `benchmarks/comparison/equivalence.json` with all observed metrics and thresholds;
- `benchmarks/comparison/execution.json` with script status and duration;
- `benchmarks/comparison/report/summary.md` for the reviewer response;
- `benchmarks/comparison/report/metrics.csv` for a supplementary table; and
- figures and biological marker-overlap tables under `benchmarks/comparison/report/artifacts`.

A failed threshold makes the final command fail but does not stop later comparisons from
running, so incomplete equivalence still yields a complete diagnostic report.

For the Theislab/HMGU GPU cluster, use the documented
[SLURM + uv + localscratch workflow](cluster/README.md).

See the complete [CPU/GPU coverage inventory](benchmarks/comparison/COVERAGE.md) for
the mapping from public methods to evidence scripts, and [`OPEN.md`](OPEN.md) for what this
evidence does **not** establish — scale limits, unreviewed thresholds, upstream issues found,
GPU portability, and what automated validation would require.

## Run time

The benchmark notebooks and source figures are available under
[`benchmarks/speed`](benchmarks/speed). Their environment and data choices should be kept
separate from the CPU/GPU equivalence thresholds above: speed measures performance, while
the comparison suite measures numerical and biological agreement.

Automatic execution requires GPU-backed CI. The manual
[`gpu-equivalence` workflow](.github/workflows/gpu-equivalence.yml) targets a self-hosted
Linux runner labeled `gpu`; otherwise submit the documented Slurm job directly.
