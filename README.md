# rapids_singlecell-reproducibility

This repository contains all scripts to reproduce the associated publication.

## Use cases

- [`use_cases/100M-working.ipynb`](use_cases/100M-working.ipynb) contains the Dask-based
  Tahoe-100M workflow.
- [`benchmarks/speed/notebooks`](benchmarks/speed/notebooks) contains matching CPU and GPU
  notebooks and source figures for the million-cell mouse-brain workflow.

## Environment

`pyproject.toml` and `uv.lock` pin the environment the benchmarks use (rapids-singlecell 0.17.0,
CUDA 13, Python 3.13) and the CPU reference packages. `uv sync` installs
it, including this repository in editable mode.

## CPU/GPU equivalence

[`benchmarks/comparison`](benchmarks/comparison) runs scanpy, squidpy, decoupler or pertpy on the
CPU and rapids-singlecell on the GPU on the same input (`scanpy_core/`, `scanpy_extended/`,
`squidpy/`, `decoupler/`, `pertpy/`, `biological_pipeline/`). Each script stores both outputs as
Zarr; [`evaluate.py`](benchmarks/comparison/evaluate.py) then computes the comparison metrics and
scores them against [`criteria.py`](benchmarks/comparison/criteria.py), the one place thresholds
are written.

```bash
uv run benchmarks/comparison/run_structured.py                               # everything
uv run benchmarks/comparison/run_structured.py scanpy_core/preprocessing.py  # a subset
EQUIVALENCE_OUT=<run> uv run benchmarks/comparison/evaluate.py             # re-score a run, no GPU
```

Everything a run writes goes to `out/`, or to `$EQUIVALENCE_OUT`. The same suite ships as a container, which needs only an
NVIDIA driver on the host:

```bash
docker build -t rsc-equivalence .
docker run --rm --gpus all -v "$PWD/out:/out" rsc-equivalence
```

### Results

[`COMPARISON_RESULTS.md`](COMPARISON_RESULTS.md) is the generated report of the last full run: hardware, package
versions, every criterion with its observed value, and every measurement recorded without a
verdict. It is the only committed output. To update it after a full run:

```bash
cp out/report/summary.md COMPARISON_RESULTS.md
```

### Criteria

- Preprocessing, HVG selection, PCA and ligrec means use `numpy.allclose` with `atol` raised to float32
  resolution; the value at numpy's default `atol` is reported alongside, without a verdict. Other
  outputs are gated on correlation or agreement.
- UMAP is held to the CPU's own seed-to-seed variation.
- Harmony compares `harmony2` with harmonypy 2, both with `max_iter_harmony=20`.
