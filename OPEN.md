# Open questions and known limitations

What the CPU/GPU equivalence work does **not** yet establish, stated plainly so a reviewer does
not have to infer it. Every number below comes from the suite itself; the generated
`report/summary.md` prints each failing criterion with the diagnosis behind it.

For context, what the suite **does** establish: on the complete 20-group run
([`snapshots/2026-07-31-expanded`](benchmarks/comparison/snapshots/2026-07-31-expanded)), 200 of 210
metrics pass, and after six targeted probes none of the 10 failures is attributable to a
rapids-singlecell defect.
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

Every threshold in the suite carries an empty `basis` field — the failing ones now carry a
`diagnosis`, which says what was measured, not that the threshold was agreed. They were all chosen
before any of this analysis and have not been reviewed with the method owners, so **"200/210
passing" reads stronger than it is**: a green metric with an arbitrary threshold is weak evidence in
the same way a red one is. The reds have at least been investigated; the greens have not.

Four criteria were diagnosed as unsatisfiable by any correct implementation and were
**deliberately left failing** rather than widened. Re-specifying them is a decision for the method
owners, not something to do to make a run green. Their supporting measurements are in the next
section.

Stochastic comparisons also report a single realization at one fixed seed. A reseeded-CPU baseline
is now recorded beside the embedding-overlap criteria as evidence, but distributions or confidence
intervals over several seeds are not yet reported.

## 3. What the four diagnoses rest on

Each failing criterion carries a one-line diagnosis in the generated report. The measurements
behind those lines:

**`normalize_total.max_abs_error` is unattainable in float32.** Observed 1.220703125e-4 against
`<= 1e-5`. pbmc3k normalized to `target_sum=10000` reaches 1751.05, where one float32 ULP is
_exactly_ 1.220703125e-4 — the "error" is the smallest difference float32 can represent there, and
the two values agree to nine digits. Relative to the data scale the disagreement is 6.9712e-8.

**`umap.*.trustworthiness` scores one implementation, not agreement.** Trustworthiness compares an
embedding to its own input. Scanpy never reaches 0.9 on pbmc3k:

| CPU seed | Trustworthiness (pbmc3k) |
| -------: | -----------------------: |
|        0 |                  0.86958 |
|        1 |                  0.86811 |
|        2 |                  0.87050 |

The GPU value, 0.86983, sits inside that range. On `pbmc68k_reduced` both criteria pass.

**`umap.cross_embedding_knn_overlap` has no reference point.** Rerunning the CPU embedding with
only the seed changed:

| Comparison               |      pbmc3k | `pbmc68k_reduced` |
| ------------------------ | ----------: | ----------------: |
| CPU seed 0 vs CPU seed 1 |     0.42800 |           0.55752 |
| CPU seed 0 vs CPU seed 2 |     0.41878 |           0.55581 |
| **CPU seed 0 vs GPU**    | **0.40190** |       **0.58848** |

On `pbmc68k_reduced` the GPU embedding is **closer** to the CPU than a reseeded CPU run is, and the
criterion still fails. The baseline is measured in-run rather than hard-coded because it is
platform-dependent: it moved from 0.5558 to 0.5859 between machines.

**`calculate_niche` measures Leiden backend choice, not correctness.** The UTAG features and their
PCA agree numerically — components 0–20 match at |r| >= 0.9998 and carry 99.997% of the variance,
and subspace alignment is >= 0.999998 for the first 15 components. What diverges is the clustering:

| On one identical graph, only the Leiden backend changed |    ARI |    NMI | Clusters |
| ------------------------------------------------------- | -----: | -----: | -------: |
| UTAG: CPU leidenalg vs CPU igraph                       | 0.5041 | 0.7120 | 12 vs 12 |
| UTAG: CPU leidenalg vs GPU cuGraph                      | 0.5941 | 0.7554 | 12 vs 14 |
| Neighborhood: CPU leidenalg vs CPU igraph               | 0.9793 |      — | 41 vs 40 |
| Neighborhood: CPU leidenalg vs GPU cuGraph              | 0.9454 |      — | 41 vs 34 |

Scanpy's own two backends agree no better with each other than the GPU agrees with Scanpy. For the
`neighborhood` flavor part of the divergence also happens _before_ Leiden, and there
rapids-singlecell is the accurate side — the profile holds only 755 distinct rows across 4668
cells, so 89.67% of cells have an exact distance tie at the k-th neighbour:

| Checked against exact float64 ground truth, first 13 non-self neighbours | Distances materially worse |  Rows | Max excess |
| ------------------------------------------------------------------------ | -------------------------: | ----: | ---------: |
| `sc.pp.neighbors` default                                                |                       2055 |   743 |     +3.014 |
| `sc.pp.neighbors` forced exact (`transformer="sklearn"`)                 |                       2055 |   743 |     +3.014 |
| `rsc.pp.neighbors` brute                                                 |                      **0** | **0** |   +4.9e-06 |

Recall against ground truth tells the same story: rsc 0.899 versus scanpy 0.858 at k=14.

### Every failure was probed for a GPU-side defect; none was found

| Suspected defect                            | Verdict                                                                                 |
| ------------------------------------------- | --------------------------------------------------------------------------------------- |
| `rsc.pp.neighbors` returns wrong neighbours | No. Exact to 4.9e-06 against float64 ground truth; the CPU reference is the inexact one |
| `n_neighbors` off-by-one in the graph       | No. Both yield 14 real neighbours for `n_neighbors=15` on clean data; Jaccard 1.000000  |
| Approximate search degrading the GPU graph  | No. GPU approximate versus GPU brute is 1.000000                                        |
| `rsc.tl.leiden` finds a worse partition     | No. Modularity 0.916633 versus leidenalg 0.917280 — 0.07% on a very flat objective      |
| `rsc.pp.pca` wrong on the UTAG features     | No. Subspace alignment >= 0.999998 over the components carrying 99.997% of the variance |
| `rsc.pp.normalize_total` numerically wrong  | No. Agrees with scanpy to exactly one float32 ULP                                       |

## 4. Upstream issues found, not yet filed

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

## 5. GPU hardware and portability

**The pinned wheel does not run on every GPU.** `rapids-singlecell-cu12==0.16.1` fails on V100
(sm_70): two runs reached CUDA and then died inside preprocessing kernels with `named symbol not
found`. A100 and H100 hosts run the same code fine. This is an architecture-coverage gap in the
build, and it is the one finding that directly undercuts reproducibility: a third party cannot
reproduce these results on arbitrary GPU hardware. The container pins the software environment; it
does not close this gap.

**A visible GPU is not always a usable one.** On one run the GPU appeared normally in `nvidia-smi`
while CUDA returned `cudaErrorNoDevice`, with an unknown MIG state and a request for a GPU reset.

[`benchmarks/comparison/gpu_smoke_check.py`](benchmarks/comparison/gpu_smoke_check.py) detects both
failure modes in seconds — it records GPU, compute capability, driver and CUDA runtime, executes
one real kernel, always writes `gpu-smoke.json`, and exits 90 rather than spending a whole run. It
distinguishes infrastructure failures from scientific ones; it does not fix either.

Consequently every number in the current snapshots comes from **A100 or H100 only**.

## 6. No automated validation

There is no GPU-backed CI, so all of this is a **point-in-time result rather than a regression
guard**. Nothing prevents a future change from silently breaking equivalence.

The manual [`gpu-equivalence` workflow](.github/workflows/gpu-equivalence.yml) targets a
self-hosted runner labelled `gpu` that does not exist yet, and has never been exercised. Making
this automatic needs, in order:

1. a self-hosted GPU runner carrying the workflow's labels;
2. one run proving that reports and per-script logs upload **even when thresholds fail**, since a
   failing suite is the normal case here and must still publish evidence;
3. an agreed blocking policy — infrastructure failures, missing result records, and a failed smoke
   check should always block, while provisional scientific thresholds should only block once
   reviewed and baselined; and
4. scheduled or release-triggered runs, once the resource cost is understood. Keep pull-request
   execution manual or scoped until then.

## 7. Data provenance is not pinned

Datasets are downloaded at run time by `scanpy.datasets` and `squidpy.datasets` with no version
pin or checksum, and the generated report does not record dataset provenance. A silent upstream
change to any of them would move the numbers without any signal.
