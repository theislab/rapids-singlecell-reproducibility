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

Each script writes a JSON record identifying the method, reference package and version,
dataset and tier, plus one entry per metric:

| Field        | Meaning                                                |
| ------------ | ------------------------------------------------------ |
| `metric`     | Metric name, unique within a method group              |
| `observed`   | Measured value                                         |
| `comparison` | `<=`, `>=`, or `observed` for a non-gating measurement |
| `tolerance`  | Threshold, or `null` for a non-gating measurement      |
| `criterion`  | Human-readable form of the threshold                   |
| `gating`     | Whether the metric decides pass/fail                   |
| `basis`      | Why the threshold is what it is; empty until reviewed  |
| `diagnosis`  | For a failing criterion, what the investigation found  |
| `passed`     | Result; always `true` for non-gating measurements      |

Metrics are split into criteria that gate the suite and measurements recorded as evidence;
`collect_results.py` rejects a record that mixes the two. **Thresholds are never widened to
make a run green.** A failing criterion keeps its value and carries its `diagnosis`, which the
generated report prints beside it. After running the comparisons individually, aggregate the
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

See the complete [CPU/GPU coverage inventory](benchmarks/comparison/COVERAGE.md) for
the mapping from public methods to evidence scripts, and [`OPEN.md`](OPEN.md) for what this
evidence does **not** establish — scale limits, unreviewed thresholds, upstream issues found,
GPU portability, and what automated validation would require.

## Container

The suite ships as a container, so reproducing it does not mean rebuilding the environment
by hand. CUDA comes from the wheels pinned in `uv.lock`; the only host requirement is an
NVIDIA driver.

```bash
docker build -t rsc-equivalence .
```

```bash
docker run --rm --gpus all -v "$PWD/out:/out" rsc-equivalence
```

Everything a run writes — result records, the report, downloaded datasets — lands in `out/`.
The image takes the same arguments as `run_structured.py`, and also carries a fast GPU check
that tells an unusable GPU apart from a scientific failure and exits `90` when
rapids-singlecell kernels cannot run at all:

```bash
docker run --rm --gpus all -v "$PWD/out:/out" rsc-equivalence /repro/benchmarks/comparison/gpu_smoke_check.py --output /out/gpu-smoke.json
```

```bash
docker run --rm --gpus all -v "$PWD/out:/out" rsc-equivalence /repro/benchmarks/comparison/run_structured.py scanpy_core/preprocessing.py
```

Rootless hosts can convert and run the same image with Apptainer:

```bash
apptainer run --nv --pwd /out -B "$PWD/out:/out" rsc-equivalence.sif
```

## Run time

The benchmark notebooks and source figures are available under
[`benchmarks/speed`](benchmarks/speed). Their environment and data choices should be kept
separate from the CPU/GPU equivalence thresholds above: speed measures performance, while
the comparison suite measures numerical and biological agreement.

Automatic execution requires GPU-backed CI. The manual
[`gpu-equivalence` workflow](.github/workflows/gpu-equivalence.yml) targets a self-hosted
Linux runner labeled `gpu`; without one, run the container on a GPU host as above.
