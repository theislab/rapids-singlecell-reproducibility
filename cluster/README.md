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
5. runs every structured comparison in a separate process, including the core Scanpy and
   end-to-end biological workflows;
6. preserves per-script logs even when thresholds fail;
7. generates a reviewer-facing Markdown report, CSV table, and figures; and
8. copies the compact evidence, report, environment, and GPU manifests back to shared storage.

Optional overrides:

```bash
sbatch --export=ALL,REPRO_SOURCE=/lustre/groups/ml01/workspace/$USER/rapids-singlecell-reproducibility cluster/run_comparisons.sbatch

# Choose a durable result destination:
sbatch --export=ALL,EQUIVALENCE_DURABLE_OUTPUT=/lustre/groups/ml01/workspace/$USER/rsc-equivalence/run-01 cluster/run_comparisons.sbatch
```

The default result directory is
`benchmarks/comparison/cluster-results/$SLURM_JOB_ID` in the shared checkout.

Each result directory contains:

- `equivalence.json`: complete machine-readable metrics and thresholds;
- `execution.json`: exit status and duration for every isolated comparison;
- `results/`: one JSON record per method group;
- `report/summary.md`: the reviewer-facing aggregate report;
- `report/metrics.csv`: a flat supplementary-table source;
- `report/artifacts/`: comparison figures and biological-pipeline tables;
- `report/logs/`: complete logs, including scripts that failed before writing a record; and
- `environment.txt` and `nvidia-smi.txt`: software and hardware provenance.

The suite intentionally continues after individual failures so every comparison contributes
evidence. The Slurm job still exits non-zero at the end when a script or threshold fails.

Automatic pull-request reporting requires a GPU-backed CI runner. A manually dispatched
GitHub Actions workflow is included for a self-hosted runner carrying the `gpu` label; the
Slurm workflow remains the fallback when the cluster is not connected to GitHub Actions.
