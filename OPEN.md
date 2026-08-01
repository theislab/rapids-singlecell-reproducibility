# Open questions and known limitations

What the CPU/GPU equivalence work does **not** yet establish, stated plainly so a reviewer does
not have to infer it. The actionable version of this list lives in
[`benchmarks/comparison/TODO.md`](benchmarks/comparison/TODO.md); the measurement evidence behind
each claim is in [`benchmarks/comparison/THRESHOLDS.md`](benchmarks/comparison/THRESHOLDS.md).

For context, what the suite **does** establish: on the complete 20-group run
([`snapshots/38938766`](benchmarks/comparison/snapshots/38938766)), 200 of 210 metrics pass, and
after six targeted probes none of the 10 failures is attributable to a rapids-singlecell defect.
The deterministic core is clean — `rank_genes_groups` 60/60, decoupler 18/18, pertpy `distance`
18/18, `hvg_pca` 19/19 — and `normalize_total` agrees with Scanpy to exactly one float32 ULP.

## 1. Validated at small scale, assumed at large scale

Every equivalence comparison runs on small data:

| Dataset                    | Size        |
| -------------------------- | ----------- |
| `scanpy.datasets.pbmc3k`   | 2,700 cells |
| `pbmc3k_processed`         | 2,638 cells |
| `pbmc68k_reduced`          | 700 cells   |
| `squidpy.datasets.imc`     | 4,668 cells |
| Synthetic pertpy/decoupler | 60–200 rows |

The workflows this repository exists to showcase — the million-cell mouse brain benchmarks in
[`benchmarks/speed`](benchmarks/speed) and the Dask-based Tahoe-100M pipeline in
[`use_cases`](use_cases) — have **no equivalence coverage at all**. Nothing tests whether the
Dask/multi-chunk code paths agree with their CPU equivalents, and chunked execution is exactly
where accumulation order and partition boundaries could diverge.

This is the largest gap. Equivalence is demonstrated at 10^3 cells and assumed at 10^6–10^8.

## 2. The thresholds are unreviewed, including the passing ones

Every threshold in the suite carries an empty `basis` field. They were chosen before any of this
analysis and have not been reviewed with the method owners, so **"200/210 passing" reads stronger
than it is**: a green metric with an arbitrary threshold is weak evidence in the same way a red one
is. The reds have at least been investigated; the greens have not.

Four criteria were diagnosed as unsatisfiable by any correct implementation and were
**deliberately left failing** rather than widened — see `THRESHOLDS.md`. Re-specifying them is a
decision for the method owners, not something to do to make a run green.

Stochastic comparisons also report a single realization at one fixed seed. A reseeded-CPU baseline
is now recorded beside the embedding-overlap criteria as evidence, but distributions or confidence
intervals over several seeds are not yet reported.

## 3. Upstream issues found, not yet filed

None of these are patched in this repository.

**scanpy — a real bug.** `sc.pp.neighbors` removes the cell itself from the k-NN result _by
position_, assuming it sorts first. With exact duplicates it does not, so a genuine nearest
neighbour is dropped and self keeps a slot: the row silently ends up with `n_neighbors - 2` real
neighbours and distances shifted one place outward. Containment is exact on
`squidpy.datasets.imc` — 1261 rows retain self, 743 are inexact against float64 ground truth, all
743 lie inside the 1261, and all 3407 rows that do not retain self are exact. Max excess distance
+3.01. CPU-only reproducer:
[`scanpy_neighbors_duplicate_bug.py`](benchmarks/comparison/squidpy/scanpy_neighbors_duplicate_bug.py)
— 240 cells, 168 of them left with 13 real neighbours instead of 14.

**rapids-singlecell — a question, not a defect.** `rsc.pp.neighbors` materializes the self-loop as
an explicit stored zero in `obsp["distances"]`, so it holds `n_neighbors` entries per row where
Scanpy holds `n_neighbors - 1`. Neighbour sets are identical once the diagonal is removed and
`obsp["connectivities"]` is unaffected. Worth noting when reporting: this convention is precisely
what makes rapids-singlecell immune to the Scanpy bug above.

**squidpy — a method-design question.** `calculate_niche(flavor="neighborhood")` clusters a feature
space holding 755 distinct rows across 4,668 cells, so 89.67% of cells have an exact distance tie
at the k-th neighbour. It also calls `sc.tl.leiden` without passing `flavor`, so the result depends
on which Leiden backend is installed — and those backends disagree substantially here (ARI 0.5041
between leidenalg and igraph on the UTAG space).

## 4. GPU hardware and portability

**The pinned wheel does not run on every GPU.** `rapids-singlecell-cu12==0.16.1` fails on the
cluster's V100 nodes: jobs `38938670` and `38938725` reached CUDA and then died inside
preprocessing kernels with `named symbol not found`. A100 (`38936727`, `38938766`, `38957063`) and
H100 (`38956547`) nodes run the same code fine. This is an architecture-coverage gap in the build,
and it is the one finding that directly undercuts reproducibility: a third party cannot reproduce
these results on arbitrary GPU hardware. It is why the checked-in Slurm job cannot yet be a
generic GPU request.

**One A100 node was itself unhealthy.** Job `38938173` saw its allocated A100 in `nvidia-smi`, but
CUDA returned `cudaErrorNoDevice`; `nvidia-smi -q` reported an unknown MIG state and asked for a
GPU reset. `gpusrv29` should be reported to the HPC team if it stays in that state.

[`cluster/gpu_smoke_check.py`](cluster/gpu_smoke_check.py) now detects both failure modes in
seconds — it records GPU, compute capability, driver and CUDA runtime, executes one real kernel,
always writes `gpu-smoke.json`, and aborts with exit 90 rather than consuming the allocation. It
distinguishes infrastructure failures from scientific ones; it does not fix either.

Consequently every number in the current snapshots comes from **A100 or H100 only**.

## 5. No automated validation

There is no GPU-backed CI, so all of this is a **point-in-time result rather than a regression
guard**. Nothing prevents a future change from silently breaking equivalence.

The manual [`gpu-equivalence` workflow](.github/workflows/gpu-equivalence.yml) targets a
self-hosted runner labelled `gpu` that does not exist yet, and has never been exercised. Making
this automatic needs, in order:

1. a self-hosted GPU runner carrying the workflow's labels, or a CI-to-Slurm bridge;
2. one run proving that reports and per-script logs upload **even when thresholds fail**, since a
   failing suite is the normal case here and must still publish evidence;
3. an agreed blocking policy — infrastructure failures, missing result records, and a failed smoke
   check should always block, while provisional scientific thresholds should only block once
   reviewed and baselined; and
4. scheduled or release-triggered runs, once the resource cost is understood. Keep pull-request
   execution manual or scoped until then.

## 6. Data provenance is not pinned

Datasets are downloaded at run time by `scanpy.datasets` and `squidpy.datasets` with no version
pin or checksum, and the generated report does not record dataset provenance. A silent upstream
change to any of them would move the numbers without any signal.
