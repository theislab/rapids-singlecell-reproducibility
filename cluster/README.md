# Running the equivalence suite on the Theislab HPC

Submit [`run_comparisons.sbatch`](run_comparisons.sbatch) from a shared-storage checkout.
The login node only submits the job; environment installation, downloads, compilation, and
comparisons all happen inside the allocated GPU job.

```bash
sbatch cluster/run_comparisons.sbatch
```

Defaults: `gpu_p` / `gpu_priority`, one available GPU, six CPUs, 90 GB RAM, and four hours. The job:

1. stages the comparison repository from shared storage to `/localscratch/$USER`;
2. puts the uv, XDG, CuPy, Numba, matplotlib, joblib, and temporary caches on local scratch;
3. bootstraps a pinned uv into local scratch when it is not already available, then creates
   the Python 3.12 environment there with `uv sync --frozen`;
4. installs the pinned prebuilt `rapids-singlecell-cu12[rapids]` wheel, avoiding a toolkit-dependent source build;
5. runs [`gpu_smoke_check.py`](gpu_smoke_check.py), which initializes CuPy and executes one
   small rapids-singlecell kernel before any comparison starts;
6. runs every structured comparison in a separate process, including the core Scanpy and
   end-to-end biological workflows;
7. preserves per-script logs even when thresholds fail;
8. generates a reviewer-facing Markdown report, CSV table, and figures; and
9. copies the compact evidence, report, environment, and GPU manifests back to shared storage.

## GPU smoke check

The smoke check separates infrastructure failures from scientific ones. It records the GPU
name, compute capability, driver version, CUDA runtime version and memory it actually sees,
then runs `rsc.pp.normalize_total` on a 64x32 matrix. Its record is always written to
`gpu-smoke.json` in the result directory, including on failure.

A failing check exits the job with code `90` before the suite starts, because a node that
cannot run rapids-singlecell kernels produces no scientific evidence and would otherwise
consume the whole allocation. Two failure modes seen on this cluster motivate it: V100 nodes
where the pinned wheel raised `named symbol not found` inside preprocessing kernels, and an
A100 node visible in `nvidia-smi` where CUDA returned `cudaErrorNoDevice`.

Optional overrides:

```bash
sbatch --export=ALL,REPRO_SOURCE=/lustre/groups/ml01/workspace/$USER/rapids-singlecell-reproducibility cluster/run_comparisons.sbatch

# Choose a durable result destination:
sbatch --export=ALL,EQUIVALENCE_DURABLE_OUTPUT=/lustre/groups/ml01/workspace/$USER/rsc-equivalence/run-01 cluster/run_comparisons.sbatch

# Rerun only selected comparisons, for iterating on one method:
sbatch --time=01:00:00 \
  --export=ALL,EQUIVALENCE_SCRIPTS='scanpy_core/preprocessing.py biological_pipeline/pbmc3k.py' \
  cluster/run_comparisons.sbatch
```

`EQUIVALENCE_SCRIPTS` accepts a space- or comma-separated subset of the paths listed in
[`run_structured.py`](../benchmarks/comparison/run_structured.py); an unknown name is an
error rather than a silent skip. A subset run keeps the untouched result records in place so
the aggregate still covers the suite, marks `execution.json` with `"partial": true`, and
prints a warning at the top of the generated report. Only a full run should be promoted to a
result snapshot.

The default result directory is
`benchmarks/comparison/cluster-results/$SLURM_JOB_ID` in the shared checkout.

Each result directory contains:

- `equivalence.json`: complete machine-readable metrics and thresholds;
- `execution.json`: exit status and duration for every isolated comparison;
- `results/`: one JSON record per method group;
- `report/summary.md`: the reviewer-facing aggregate report;
- `report/metrics.csv`: a flat supplementary-table source;
- `report/artifacts/`: comparison figures and biological-pipeline tables;
- `report/logs/`: complete logs, including scripts that failed before writing a record;
- `gpu-smoke.json`: the post-install GPU check, including its failure details; and
- `environment.txt` and `nvidia-smi.txt`: software and hardware provenance.

The suite intentionally continues after individual failures so every comparison contributes
evidence. The Slurm job still exits non-zero at the end when a script or threshold fails.

Automatic pull-request reporting requires a GPU-backed CI runner. A manually dispatched
GitHub Actions workflow is included for a self-hosted runner carrying the `gpu` label; the
Slurm workflow remains the fallback when the cluster is not connected to GitHub Actions.
