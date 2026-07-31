# Reproducibility follow-ups

The comparison harness is intentionally allowed to record failing equivalence
thresholds. A red metric is useful evidence; a script that crashes before writing
its structured result is not.

## Immediate

- [ ] Collect Slurm job `38938766` after it finishes and inspect
  `execution.json`, `equivalence.json`, `report/summary.md`, and every failed
  script log. Fix only missing-result or runtime/API failures first; do not relax
  numerical thresholds merely to make the run green.
- [ ] Promote a complete expanded-suite run from `cluster-results/<job-id>` into
  a versioned result snapshot after all 20 scripts produce records. Keep the
  existing `temp-results/38936727` snapshot as provenance for the earlier
  15-script suite.
- [ ] Revisit `calculate_niche` last. Its earlier comparison produced a valid
  record but failed 6 of 9 metrics, so it needs method-specific investigation
  rather than a broad tolerance increase.

## Cluster portability

- [ ] Add a fast post-install smoke check that initializes CuPy and runs one
  small rapids-singlecell kernel before starting the full suite. Persist the GPU,
  driver, CUDA, compute-capability, and failure details when this check fails.
- [ ] Resolve compatibility of the locked `rapids-singlecell-cu12==0.16.1` wheel
  with the cluster's V100 nodes. Jobs `38938670` and `38938725` reached CUDA but
  RSC preprocessing kernels failed with `named symbol not found`; the A100 node
  used by jobs `38936727` and `38938766` runs the kernels. Prefer a wheel with the
  required architecture coverage so the checked-in job can remain a generic GPU
  request.
- [ ] Report `gpusrv29` to the HPC team if it remains unhealthy. Job `38938173`
  saw the allocated A100 in `nvidia-smi`, but CUDA returned `cudaErrorNoDevice`;
  `nvidia-smi -q` reported an unknown MIG state and requested a GPU reset.

## Scientific validation

- [ ] Review every tolerance with method owners and record its scientific basis.
  Keep deterministic numerical tolerances separate from stochastic structural
  metrics such as graph overlap, ARI/NMI, and embedding trustworthiness.
- [ ] Run stochastic comparisons over multiple fixed seeds and report
  distributions or confidence intervals, not only a single realization.
- [ ] Add at least one larger real dataset after the PBMC3k biological pipeline
  is stable, then report label-transfer accuracy, marker recovery, clustering
  agreement, and runtime/memory alongside low-level numerical metrics.
- [ ] Pin or checksum every externally downloaded dataset and record dataset
  provenance in the generated summary.
- [ ] Have an independent maintainer review the CPU reference choice and matched
  parameters for each method, especially Harmony and downstream biological
  metrics.

## Automation and reviewer evidence

- [ ] Configure a self-hosted GPU runner with the labels used by
  `.github/workflows/gpu-equivalence.yml`, then exercise the manual workflow and
  verify that reports and per-script logs upload even when thresholds fail.
- [ ] Add scheduled or release-triggered GPU runs once runner capacity and data
  access are reliable. Keep pull-request execution manual or scoped until the
  resource cost is understood.
- [ ] Link a stable generated `report/summary.md`, metric CSV, environment lock,
  GPU inventory, and raw logs from the reviewer response. State explicitly that
  automatic enforcement requires GPU CI.
- [ ] Decide which failures block releases: infrastructure and missing records
  should always block; provisional scientific thresholds should become blocking
  only after they are reviewed and baselined.

## Current run handoff

Job `38938766` is running on `gpusrv26` from the `gpu-equivalence` branch. Its
durable evidence will be written to:

```text
/lustre/groups/ml01/workspace/selman.ozleyen/rapids-singlecell-reproducibility/benchmarks/comparison/cluster-results/38938766
```

Check it with:

```bash
ssh hpc-submit01.scidom.de \
  'sacct -j 38938766 --format=JobID,NodeList,State,Elapsed,ExitCode'
```
