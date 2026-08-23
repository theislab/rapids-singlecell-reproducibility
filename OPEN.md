# Open questions and known limitations

What the CPU/GPU equivalence work does **not** yet establish, stated plainly so a reviewer does not
have to infer it.

**Where the numbers live.** Counts and per-metric values belong to a run, so this document does not
restate them. The current evidence is [`EVIDENCE.md`](EVIDENCE.md), which lists every gating
criterion and every recorded measurement with its observed value, plus the measured diagnosis
behind each failure, and dates its own measurements. It is the only committed record of a run —
the run directory itself is not kept, so the numbers there are checked by reading them, and
reproduced by rerunning the container rather than by re-scoring a stored copy:

```bash
docker run --rm --gpus all -v "$PWD/out:/out" rsc-equivalence
```

What follows is the interpretation: what the suite covers, what it cannot show, and why each red
metric is red. Only quantities that do not change from run to run are quoted here.

For context, what the suite **does** establish: on that run **11 criteria fail out of 150**, and none
of the failures is attributable to a rapids-singlecell defect. They come from three causes: three `allclose`
criteria that float32 cannot satisfy at all, one that is a genuine relative difference too small to
move any downstream result (both in section 3), and seven stochastic criteria that ask for more
agreement than the CPU reference shows against itself (section 4).

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

## 2. Most thresholds are still unreviewed

The deterministic criteria now state the standard the manuscript declares — `numpy.allclose` at
default parameters (`rtol=1e-5`, `atol=1e-8`) — reported as `allclose_excess`, the worst
elementwise difference as a fraction of that envelope, where `<= 1` means the two arrays are
`allclose`. Harmony's per-component correlation floor of 0.95 is likewise the manuscript's own
Harmony criterion. Those criteria carry a `basis` naming the published standard.

That standard gates only where the publication actually asserts it. The Methods name
normalization, HVG selection and PCA; the abstract adds log transformation and scaling. Nothing
else is claimed, so nothing else is judged: spatial statistics, activity inference, perturbation
signatures, pertpy distances, `regress_out`, `score_genes` and the QC metrics are measured
against the same envelope and reported without a verdict, because a red row traceable to no
published claim is not evidence of anything. Those 31 comparisons are still computed from the
stored arrays — `criteria.py` lists them as `EVIDENCE` — and five of them sit outside the
envelope, assessed in [`NUMERICAL_VALIDATION.md`](NUMERICAL_VALIDATION.md) alongside the four
that gate. Storing the arrays is what makes this cheap: demoting a criterion costs no
information, because the number is still there to be argued about later.

**No other threshold does.** The correlation floors, the Jaccard floors, and every ARI/NMI
threshold were chosen before any of this analysis and have not been reviewed with the method
owners. For stochastic methods the manuscript names adjusted Rand index and "preservation of local
neighborhood structure" but states **no threshold at all**, so those numbers have no external
source. "N of M passing" therefore reads stronger than it is: a green metric with an arbitrary
threshold is weak evidence in the same way a red one is.

The suite also gates only on _differences_. Absolute quality scores of a single implementation —
UMAP trustworthiness, per-backend annotation accuracy, per-backend cell-type NMI — are recorded as
evidence rather than asserted, because they do not test CPU/GPU agreement at all and the CPU
reference can fail them on its own.

Stochastic comparisons still report a single realization at one fixed seed. A reseeded-CPU baseline
is recorded beside the embedding-overlap criteria as evidence, but distributions or confidence
intervals over several seeds are not yet reported.

## 3. The declared validation standard is not met by every operation it names

The manuscript Methods state that deterministic operations agree within `numpy.allclose` at default
parameters. The suite states that criterion directly, as `allclose_excess`: the worst elementwise
difference divided by `numpy.allclose`'s own envelope, `atol + rtol * |b|`, so `<= 1` means the two
arrays are `allclose`. Of the 50 comparisons measured this way, nine exceed the envelope — four of
them under a criterion that gates, five on operations the publication makes no claim about. Each
records the reference magnitude at its worst element, so which of the criterion's two terms decided
the verdict is measured rather than inferred.

Some are decided by `atol=1e-8`, an absolute floor calibrated for float64, at elements where the
quantity passes through zero — **no float32 implementation can satisfy it there**, whatever the
implementation. The rest are decided by the relative term and are genuine disagreements above
`rtol=1e-5`, small enough to move no downstream result but larger than the tolerance the publication
cites. None is a rapids-singlecell defect: where a paired Pearson correlation is measured on the same
arrays, it is perfect.

**[`NUMERICAL_VALIDATION.md`](NUMERICAL_VALIDATION.md) is the assessment**, with the per-comparison
table, the split between the two causes, and the cross-architecture drift. Deciding what the Methods
should say is an authors' decision, not something to settle by adjusting a threshold here.

## 4. What the remaining red criteria rest on

Seven stochastic criteria fail: the cross-embedding overlap on both UMAP comparisons, and five
`calculate_niche` criteria. Every one is still in place at its original value.

Two criteria that used to fail here no longer gate. `umap.cpu.trustworthiness` and
`umap.gpu.trustworthiness` score a single embedding against its own input, so they never tested
CPU/GPU agreement; they are recorded as evidence instead. That was a deliberate re-specification, not
a relaxation to reach green — and it does leave a gap, noted at the end of this section.

**`umap.cross_embedding_knn_overlap` has no reference point.** UMAP is stochastic, so a CPU-vs-GPU
overlap is only interpretable next to how far the CPU reference is from itself. The suite therefore
reruns the CPU embedding with only the seed changed and records that baseline beside the criterion,
as `umap.cpu_reseeded_knn_overlap` and `umap.cross_embedding_overlap_vs_cpu_baseline` in
[`EVIDENCE.md`](EVIDENCE.md).

On `pbmc68k_reduced` the GPU embedding has been **closer** to the CPU embedding than a reseeded CPU
run is, while the criterion still fails. The baseline is measured in-run rather than hard-coded
precisely because it is platform-dependent and moves between machines — which is also why no value
for it is quoted here.

**`calculate_niche` measures Leiden backend choice, not correctness.** The UTAG features and their PCA
agree numerically — the leading components match at near-unit correlation, carry essentially all of
the variance, and their subspace alignment is near-exact. What diverges is the clustering.

Hold the graph identical and change only the Leiden backend, and **Scanpy's own two backends
(leidenalg and igraph) agree no better with each other than the GPU agrees with Scanpy** — on the
UTAG feature space they agree well below the threshold this criterion asserts. So the criterion is
measuring backend choice and tie-breaking, not correctness.

For the `neighborhood` flavor part of the divergence happens _before_ Leiden, and there
**rapids-singlecell is the accurate side**: checked against exact float64 ground truth, the GPU
neighbour graph is exact while `sc.pp.neighbors` returns strictly worse neighbours on a substantial
minority of rows, whether or not exact search is forced. The cause is that the neighborhood profile
holds far fewer distinct rows than cells, so the large majority of cells have an exact distance tie
at the k-th neighbour and the partition is not well determined.

Both findings are reproduced by
[`niche_divergence_diagnostic.py`](benchmarks/comparison/squidpy/niche_divergence_diagnostic.py),
which separates features, kNN graph and Leiden backend; the figures are its output rather than the
suite's, and are deliberately not copied here.

### The gap left by demoting the quality scores

Nothing now gates UMAP embedding quality, per-backend annotation accuracy, or per-backend cell-type
NMI. Each was an absolute score of one implementation, so none of them tested equivalence — but their
removal means an upstream regression that degraded **both** backends equally would leave
`accuracy_difference` near zero and the suite green while the pipeline produced nonsense. The paired
difference criteria still gate for accuracy and cell-type NMI; for trustworthiness, the paired
difference is itself recorded rather than gating, so that one has no backstop at all.

### Every failure was probed for a GPU-side defect; none was found

| Suspected defect                            | Verdict                                                                     |
| ------------------------------------------- | --------------------------------------------------------------------------- |
| `rsc.pp.neighbors` returns wrong neighbours | No. Exact against float64 ground truth; the CPU reference is the inexact one |
| `n_neighbors` off-by-one in the graph       | No. Both yield the same real neighbours on clean data; Jaccard is exact      |
| Approximate search degrading the GPU graph  | No. GPU approximate versus GPU brute agree exactly                          |
| `rsc.tl.leiden` finds a worse partition     | No. Modularity differs by a fraction of a percent on a very flat objective   |
| `rsc.pp.pca` wrong on the UTAG features     | No. Near-exact subspace alignment over the components carrying the variance  |
| `rsc.pp.normalize_total` numerically wrong  | No. Agrees with scanpy to one float32 ULP                                    |

## 5. Upstream issues found, not yet filed

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

## 6. GPU hardware and portability

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

Consequently every number in the reported evidence comes from **A100 or H100 only**.

## 7. No automated validation

Nothing runs the suite on a schedule, on a pull request, or on a release. Every result in this
repository was produced by a person running the container on a GPU host, so all of it is a
**point-in-time result rather than a regression guard**. Nothing prevents a future change from
silently breaking equivalence.

What that costs is not hypothetical. The original one-method scripts are not run by anything, and
`umap/umap.py` accumulated two independent defects that survived undetected until 2026-08-23: the
filename shadowed the installed `umap` package, so `import umap.umap_` inside umap-learn failed
before any comparison ran, and behind that a typo (`rust_sc` assigned, `trust_sc` asserted) raised
`NameError`. The script had therefore never executed successfully at all. Both are fixed, and the
file is now `umap/umap_embedding.py`; the point is that nothing but a manual run would have found
either.

Making it automatic needs, in order:

1. a GPU runner the project actually controls, since the suite cannot run on hosted CI;
2. one run proving that reports and per-script logs are published **even when criteria fail**, since
   a failing suite is the normal case here and must still produce evidence;
3. an agreed blocking policy — infrastructure failures, missing result records, and a failed smoke
   check should always block, while provisional scientific thresholds should only block once
   reviewed and baselined; and
4. scheduled or release-triggered runs, once the resource cost is understood.

Evaluation is cheaper than measurement and needs no GPU, so re-scoring a run's stored records against
changed criteria is the one part of this that could be automated today.

## 8. One Methods statement is out of date

The Methods say `rapids_singlecell.ptg.Distance` "currently supports the E-distance metric".
Version 0.16.1 — the version pinned here and benchmarked in the manuscript — declares nine in
`SUPPORTED_METRICS`, each with its own GPU implementation: `edistance` through a CUDA kernel,
`wasserstein` through Sinkhorn, and seven pseudobulk metrics computed from group mean vectors.
All nine are compared against pertpy here and agree. The sentence understates the package.

The same paragraph describes `bootstrap`, `onesided_distances` and multi-GPU aggregation.
`onesided_distances` is now covered; `bootstrap` and the multi-GPU path are not, and
`create_contrasts`/`validate_contrasts`/`contrast_distances` are not mentioned in the Methods at
all. Section 10 below lists the gaps.

## 9. Data provenance is not pinned

The six `scanpy.datasets` and `squidpy.datasets` loaders fetch whatever upstream serves at run
time, with no version pin or checksum, and the generated report does not record dataset
provenance. A silent upstream change to any of them would move the numbers without any signal.

One exception, so this is not read more broadly than it holds: the Harmony comparison retrieves its
two inputs through `pooch` with `known_hash` md5 pins, so that comparison would fail loudly rather
than silently if its data changed.

## 10. API surface that is not compared at all

`EVIDENCE.md` names every method group that *is* compared, so this lists only the gaps. They were
tracked in a separate coverage inventory, which is gone — one generated file plus this document is
the whole of the written record now, and a hand-maintained table of what the suite covers went stale
against it immediately.

- **`ptg.Distance.bootstrap`** — output is stochastic and the manuscript states no agreement
  criterion for it, so there is nothing to gate against.
- **`ptg.Distance` contrast API** — `create_contrasts`, `validate_contrasts`, `contrast_distances`.
- **Multi-GPU execution.** `multi_gpu` appears exactly twice in the suite, both in
  `pertpy/distance.py` and both passed `False`; no other comparison exercises the parameter at all.
  The device-splitting and cross-device aggregation paths described in the Methods are therefore
  **untested**. This is the largest single gap here: the code path a multi-GPU user takes has no
  evidence behind it.
- **Squidpy's `spatialleiden` flavor**, which rapids-singlecell does not implement.

`kmeans` *is* compared, and is worth noting because it is exported by `rapids_singlecell.tl` without
appearing in `docs/api/scanpy_gpu.md`.

## 11. The original one-method scripts were removed

The repository previously carried 18 one-method scripts (`hvg/`, `normalize/`, `pca/`, `umap/`,
`leiden/`, and so on), each running one method on CPU and GPU and asserting a single
`numpy.assert_allclose`. They are gone, and this records why rather than leaving it to the diff.

Every method they covered has a gating criterion in `EVIDENCE.md`, so no coverage was lost. The
reason for removing rather than keeping them is that their tolerances were hand-picked and
**looser than the standard the Methods declare**: `scale` and `normalize_pearson_residuals`
asserted at `atol=1e-6`, `regress_out` at `atol=1e-5`, against declared `numpy.allclose` defaults
of `rtol=1e-5, atol=1e-8`. The suite's failures for those same operations are decided by the
absolute term, so those scripts passed only because their absolute floor was 100x to 1000x wider
than the one the manuscript states. Keeping them would have left two verdicts on the same
operations in the same repository, the greener one measured against a bar that appears nowhere in
the paper. The measured excesses are in [`EVIDENCE.md`](EVIDENCE.md) and are not restated here.

Their history is intact in git; they were contributed in the two merged benchmark pull requests.

## 12. Harmony has two algorithm flavors, and the default is not the harmonypy one

`rsc.pp.harmony_integrate` takes `flavor: Literal["harmony2", "harmony1"] = "harmony2"`. Harmonypy
implements the harmony1 algorithm, so **the rapids-singlecell default is not the CPU reference's
algorithm** — comparing them without pinning `flavor` compares two different methods.

Measured on the harmonypy 3,500-cell donor benchmark, worst relative L2 per component:

| Comparison | Worst relative L2 |
| --- | --- |
| harmonypy 0.2.0 run now, against the repository's stored harmonized file | 0.053 |
| rapids-singlecell `flavor="harmony1"`, against harmonypy run now | 0.053 |
| rapids-singlecell at its **default** `flavor="harmony2"`, against the stored file | **0.523** |

The first row matters: the stored reference is sound, reproduced by harmonypy today to the same
tolerance the GPU achieves. So the tenfold divergence in the last row is the flavor, not drift in
the reference and not a GPU defect. The suite pins `flavor="harmony1"` and agrees; the original
one-method script did not pin it and disagreed.

Two consequences worth stating. Any Harmony comparison must pin `flavor`, or it silently measures
algorithm choice — the same failure mode as `sc.tl.leiden` without `flavor` in section 13 and in
`calculate_niche`. And the published criterion is a per-component Pearson correlation, which is
invariant to a per-component scale factor; relative L2 is not, so it is gated alongside at the
threshold the original script used, to keep a magnitude difference from passing unseen.

## 13. The one-method scripts disagree with the suite about Leiden

`leiden/leiden.py` and the suite gate the same metric at the same threshold — Leiden ARI >= 0.9 —
and reach opposite verdicts. The suite passes on `pbmc68k_reduced`; the standalone script fails on
pbmc3k at ARI 0.874, with NMI 0.906 passing and both sides finding 7 clusters.

The cause is the same one already documented for `calculate_niche` in section 4: `sc.tl.leiden`
is called without `flavor`, so the CPU side runs leidenalg while the GPU side runs cuGraph. The
comparison is between two different community-detection implementations, not between CPU and GPU
executions of one. ARI 0.874 with an identical cluster count is what backend disagreement looks
like, not what a numerical defect looks like.

The threshold has been left at 0.9. It is the dataset and the unset `flavor` that make the two
disagree, and choosing which of them is the reference is a specification decision, not a
tolerance to adjust.

