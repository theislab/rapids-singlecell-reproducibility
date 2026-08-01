# Diagnosis of the failing equivalence criteria

**No criterion in this suite is ever widened, relaxed, or replaced to make a run green.** A red
metric is trustworthy evidence, and deciding that a criterion is genuinely wrong belongs to the
maintainer and the method owners. This document records what was measured and what the failures
mean; every original threshold is still in place and still failing.

Additional quantities measured while diagnosing are recorded **next to** the failing criteria as
non-gating measurements, never in place of them.

## Metric record schema

| Field        | Meaning                                                |
| ------------ | ------------------------------------------------------ |
| `metric`     | Metric name, unique within a method group              |
| `observed`   | Measured value                                         |
| `comparison` | `<=`, `>=`, or `observed` for a non-gating measurement |
| `tolerance`  | Threshold, or `null` for a non-gating measurement      |
| `criterion`  | Human-readable form of the threshold                   |
| `gating`     | Whether the metric decides pass/fail                   |
| `basis`      | What the measurement is for                            |
| `passed`     | Result; always `true` for non-gating measurements      |

`collect_results.py` rejects records that mix these — a non-gating metric carrying a tolerance,
or reporting a failure, is a malformed record. Non-gating entries exist to carry evidence, and
must never be used to demote a criterion that fails.

## Runs behind this document

| Job        | Node | What it produced                                                      |
| ---------- | ---- | --------------------------------------------------------------------- |
| `38938766` | A100 | Complete 20-group suite: 16/20 groups, 200/210 metrics, 20/20 records |
| `38956547` | H100 | Reran the three affected groups and measured the evidence quantities  |
| `38956728` | H100 | `calculate_niche` decomposition, round one                            |
| `38956767` | H100 | `calculate_niche` follow-up: PCA degeneracy and kNN graph             |
| `38956823` | A100 | Reproduced the `calculate_niche` findings on a second node            |
| `38957063` | A100 | Confirmed the reverted criteria fail again with the original values   |
| `38957332` | A100 | kNN recall against exact float64 ground truth; Leiden modularity      |
| `38957359` | A100 | Self-loop convention and correctly aligned distance comparison        |
| `38957415` | A100 | Forced scanpy exact; counted duplicate rows in the feature space      |

Job `38956547` ran while three criteria were temporarily replaced; those replacements have been
reverted, so its **PASS verdict does not apply** and is not cited as one. Its measurements of the
evidence quantities are valid and are used below.

## `normalize_total.max_abs_error`: unattainable in float32

Observed 1.220703125e-4 against `<= 1e-5`. **Still failing.**

pbmc3k normalized to `target_sum=10000` reaches a maximum value of 1751.05. One float32 ULP at
that magnitude is **exactly 1.220703125e-4** — the observed "error" is the smallest difference
float32 can represent there, and the two values agree to nine digits. Expressed relative to the
data scale the disagreement is 6.9712e-8, below one float32 epsilon.

No float32 implementation, CPU or GPU, can meet a 1e-5 absolute tolerance on data of this
magnitude. `normalize_total.float32_ulp_at_max` and `normalize_total.max_rel_error` are recorded
beside the criterion so a reviewer can see this without recomputing it. Whether the criterion
should be scale-relative is a decision for the method owners.

## `umap.*.trustworthiness`: measures one implementation, not agreement

Observed CPU 0.8699 and GPU 0.8698 against `>= 0.9`. **Still failing, including on the CPU
reference.**

Trustworthiness scores an embedding against its own input, so it does not test CPU/GPU
agreement at all. A CPU-only calibration of the biological pipeline over three UMAP seeds:

| Seed | CPU trustworthiness (pbmc3k) |
| ---: | ---------------------------: |
|    0 |                      0.86958 |
|    1 |                      0.86811 |
|    2 |                      0.87050 |

Scanpy never reaches 0.9 on this dataset. The GPU value, 0.86983, sits inside the CPU's own
seed-to-seed range. `umap.trustworthiness_difference` is recorded alongside: 4.5e-5 on the A100
node, 0.0018 on pbmc3k and 0.0032 on `pbmc68k_reduced` on the H100 node — the same order as the
CPU's own spread of ~0.0024.

On `pbmc68k_reduced` these two criteria pass (0.919 and 0.923); only the pbmc3k pipeline fails
them.

## `umap.cross_embedding_knn_overlap`: no reference point

Observed 0.4019 against `>= 0.6` on pbmc3k, and 0.5885 against `>= 0.65` on `pbmc68k_reduced`.
**Both still failing.**

UMAP is stochastic, so a CPU-vs-GPU overlap is only interpretable next to how far the CPU
reference is from itself. Rerunning the CPU embedding with only the seed changed:

| Comparison               |      pbmc3k | `pbmc68k_reduced` |
| ------------------------ | ----------: | ----------------: |
| CPU seed 0 vs CPU seed 1 |     0.42800 |           0.55752 |
| CPU seed 0 vs CPU seed 2 |     0.41878 |           0.55581 |
| **CPU seed 0 vs GPU**    | **0.40190** |       **0.58848** |

On `pbmc68k_reduced` the GPU embedding is **closer** to the CPU embedding than a reseeded CPU run
is, yet the criterion fails. On pbmc3k the GPU is 0.017 further than the weaker CPU baseline.

`umap.cpu_reseeded_knn_overlap` and `umap.cross_embedding_overlap_vs_cpu_baseline` are recorded
beside the criterion. The baseline is measured in-run rather than hard-coded because it is
platform-dependent: it moved from 0.5558 on the calibration machine to 0.5859 on the H100 node.

## `calculate_niche`: divergence located, nothing changed

`squidpy/calculate_niche.py` fails 5 of 9 metrics.
[`squidpy/niche_divergence_diagnostic.py`](squidpy/niche_divergence_diagnostic.py) separates the
features, the kNN graph, and the Leiden backend. `cellcharter` passing at ARI 0.959 was the first
clue: it ends in a Gaussian mixture, while `neighborhood` and `utag` both end in Leiden — and
Squidpy calls `sc.tl.leiden` without a flavor, so the CPU path uses leidenalg while the GPU path
uses cuGraph.

### The features are equivalent

The UTAG feature matrix and its PCA agree numerically. Components 0–20 match at |r| >= 0.9998 and
carry 99.997% of the variance; the 11 components with |r| < 0.99 each carry 0.0001–0.0006%, where
near-tied eigenvalues leave individual directions unidentified. Subspace alignment — the
identifiable quantity — is >= 0.999998 for the first 15 components on both nodes. Nothing
numerical is wrong upstream of clustering.

### UTAG: the criterion exceeds the reference's self-consistency

On one identical graph, with only the Leiden backend changed:

| Comparison                                          |        ARI |        NMI | Clusters |
| --------------------------------------------------- | ---------: | ---------: | -------: |
| CPU leidenalg vs CPU igraph                         |     0.5041 |     0.7120 | 12 vs 12 |
| CPU leidenalg vs GPU cuGraph                        |     0.5941 |     0.7554 | 12 vs 14 |
| **End-to-end CPU vs GPU (the suite's measurement)** | **0.5037** | **0.7208** | 12 vs 12 |

Scanpy's own two Leiden backends agree no better with each other (0.5041) than the GPU agrees
with Scanpy (0.5037). The UTAG feature space has no stable Leiden partition at this resolution.
The CPU-vs-CPU figure reproduces to six digits across nodes, while the end-to-end CPU-vs-GPU
value moved between 0.5037, 0.5234 and 0.5546 over three runs.

### Neighborhood: two mechanisms, neither numerical

| Comparison                                    |    ARI | Clusters |
| --------------------------------------------- | -----: | -------: |
| CPU leidenalg vs CPU igraph, identical graph  | 0.9793 | 41 vs 40 |
| CPU leidenalg vs GPU cuGraph, identical graph | 0.9454 | 41 vs 34 |
| End-to-end CPU vs GPU                         | 0.8477 | 41 vs 50 |

The end-to-end result is worse than the same-graph result, so part of the divergence happens
before Leiden. On the identical 11-dimensional scaled profile the kNN graphs differ: distances
Jaccard 0.8212, connectivities 0.8146, and it is not approximate search on either side.

**The cause is the feature space, and rapids-singlecell is the accurate side.** The profile holds
only 755 distinct rows across 4668 cells — 3913 duplicates — so 89.67% of cells have an _exact_
distance tie at the k-th neighbour boundary, against 0.00% on `pbmc68k_reduced` PCA where the two
implementations agree perfectly (Jaccard 1.000000). Checked against exact float64 ground truth on
the first 13 non-self neighbours per row (job `38957415`):

| Implementation                                           | Neighbour distances materially worse than exact | Rows affected | Max excess |
| -------------------------------------------------------- | ----------------------------------------------: | ------------: | ---------: |
| `sc.pp.neighbors` default                                |                                            2055 |           743 |     +3.014 |
| `sc.pp.neighbors` forced exact (`transformer="sklearn"`) |                                            2055 |           743 |     +3.014 |
| `rsc.pp.neighbors` brute                                 |                                           **0** |         **0** |   +4.9e-06 |

rsc matches exact ground truth to 4.9e-06; the CPU reference returns strictly worse neighbours on
743 of 4668 rows, and forcing scanpy exact reproduces its default result byte for byte
(Jaccard 1.000000), so an approximation setting is not the cause. Recall against ground truth tells
the same story: rsc 0.899 versus scanpy 0.858 at k=14.

The mechanism is confirmed (job `38957612`). `sc.pp.neighbors` stores `n_neighbors - 1 = 14`
entries per row and removes the cell itself **by position**, which only works if self sorts first
in the k-NN result. With exact duplicates it frequently does not, so a genuine nearest neighbour is
dropped instead and self keeps a slot. Affected rows are silently left with 13 real neighbours and
a distance vector shifted one place outward.

| Rows on the imc profile              |       Count |
| ------------------------------------ | ----------: |
| Retain self in `obsp["distances"]`   |        1261 |
| Inexact against float64 ground truth |         743 |
| **Inexact _and_ retain self**        |     **743** |
| **Inexact _without_ retaining self** |       **0** |
| Do not retain self, and are exact    | 3407 / 3407 |

Containment is exact: every inexact row retains self, and every row that does not retain self is
correct. The remaining 518 self-retaining rows lose a slot but sit inside a block of
zero-distance duplicates, so losing it costs them nothing measurable. `pbmc68k_reduced`, with no
duplicates, has 0 self-retaining and 0 inexact rows.

rapids-singlecell stores `n_neighbors` entries, so even on the rows where it also keeps self (2016
of them) it still carries 14 real neighbours — the convention difference noted below is what makes
it robust here. **This is a scanpy defect, not an rsc one**, and it reproduces without any GPU:
240 cells built from 12 distinct rows repeated 20 times leaves 168 of 240 rows with 13 real
neighbours instead of 14.

So the `neighborhood` criterion is failing because the suite compares an exact GPU result against
a CPU reference that is itself inexact on heavily duplicated input. This is **not** a
rapids-singlecell defect, and the earlier reading of it as one was wrong.

### Cluster-count difference

`cluster_count_difference <= 1` cannot hold across Leiden implementations. On identical input
cuGraph produced 34 clusters against leidenalg's 41, and Scanpy's own two backends differ by 1 —
already at the limit. `resolution` does not carry the same meaning across implementations.

## What was looked for on the rapids-singlecell side, and what was found

Every failure above was probed for an actual GPU-side defect (jobs `38957332`, `38957359`,
`38957415`). **No numerical or algorithmic defect was found.**

| Suspected defect                            | Verdict                                                                                 |
| ------------------------------------------- | --------------------------------------------------------------------------------------- |
| `rsc.pp.neighbors` returns wrong neighbours | No. Exact to 4.9e-06 against float64 ground truth; the CPU reference is the inexact one |
| `n_neighbors` off-by-one in the graph       | No. Both yield 14 real neighbours for `n_neighbors=15` on clean data; Jaccard 1.000000  |
| Approximate search degrading the GPU graph  | No. GPU approximate versus GPU brute is 1.000000                                        |
| `rsc.tl.leiden` finds a worse partition     | No. Modularity 0.916633 versus leidenalg 0.917280 — 0.07% on a very flat objective      |
| `rsc.pp.pca` wrong on the UTAG features     | No. Subspace alignment >= 0.999998 over the components carrying 99.997% of the variance |
| `rsc.pp.normalize_total` numerically wrong  | No. Agrees with scanpy to exactly one float32 ULP                                       |

One genuine difference, low severity and worth reporting upstream rather than working around:

**`obsp["distances"]` carries an explicit stored zero on the diagonal.** On `pbmc68k_reduced`,
`rsc.pp.neighbors` stores a diagonal entry for 700/700 rows while `sc.pp.neighbors` stores none, so
rsc's `distances` holds `n_neighbors` entries per row where scanpy holds `n_neighbors - 1` (nnz
10500 versus 9800, ratio exactly 15/14). The real neighbour sets are identical once the diagonal is
removed, and it does **not** propagate to `obsp["connectivities"]` — both are 15900 nnz with no
diagonal — so Leiden and UMAP are unaffected. The risk is downstream code that counts nonzeros per
row or assumes scanpy's convention.

## Still to review, with the method owners

- Every threshold in the suite carries an empty `basis`. The failures above are diagnosed but
  none of the criteria have been re-specified, which is deliberate.
- Whether deterministic tolerances should be scale-relative rather than absolute.
- Whether stochastic criteria should be expressed against a measured baseline. The reseeded-CPU
  overlap is recorded as evidence; making it a criterion would be a specification change.
- Stochastic comparisons still report a single realization; reporting distributions over several
  fixed seeds remains outstanding.
