# Reproducibility follow-ups

The comparison harness is intentionally allowed to record failing equivalence
thresholds. A red metric is useful evidence; a script that crashes before writing
its structured result is not.

**Thresholds are never widened, relaxed, or replaced to make a run green, and defects on the
rapids-singlecell side are not patched here.** Failures get diagnosed, recorded, and raised with
the method owners; re-specifying a criterion is their decision. See
[`THRESHOLDS.md`](THRESHOLDS.md).

## Done

- [x] Collected Slurm job `38938766`. All 20 scripts produced result records and the
      four non-zero exits were threshold assertions raised after the record was written:
      no runtime failures, no API failures, no missing records.
- [x] Promoted it to [`snapshots/38938766`](snapshots/38938766). `temp-results/38936727`
      is kept as provenance for the earlier 15-script suite and now points forward.
- [x] Diagnosed all three non-`calculate_niche` failures and **left every criterion in place**,
      writing the evidence up in [`THRESHOLDS.md`](THRESHOLDS.md): the `normalize_total` absolute
      tolerance sits below one float32 ULP at the data's magnitude; the two UMAP trustworthiness
      criteria score a single implementation rather than agreement, and the CPU reference fails
      them too; the embedding-overlap criteria demand more agreement than the CPU shows against
      itself under reseeding.
- [x] Split metric records into gating criteria and recorded measurements. The evidence
      quantities measured while diagnosing are recorded beside the failing criteria, never in
      place of them; `collect_results.py` rejects a non-gating metric that carries a tolerance
      or reports a failure.
- [x] Added [`cluster/gpu_smoke_check.py`](../../cluster/gpu_smoke_check.py), run after
      installation and before the suite. It records GPU name, compute capability, driver,
      CUDA runtime and memory, runs one rapids-singlecell kernel, always writes
      `gpu-smoke.json`, and aborts the job with exit 90 rather than burning the allocation.
- [x] Added subset selection to `run_structured.py` (`EQUIVALENCE_SCRIPTS` or positional
      arguments) so one method can be iterated without a full-suite run. Subset runs mark
      `execution.json` as partial and say so at the top of the report.
- [x] Located the source of the `calculate_niche` divergence with
      [`squidpy/niche_divergence_diagnostic.py`](squidpy/niche_divergence_diagnostic.py), which
      separates features, kNN graph, and Leiden backend. Conclusions are below; the thresholds
      were deliberately left unchanged.

## Immediate

- [ ] Take one full-suite run now that the evidence measurements have been added, and promote it
      as the current snapshot. The suite is still expected to fail: 10 metrics across four groups
      remain red by design. `snapshots/38938766` is the only complete run so far.
- [ ] Bring the four diagnosed failures to the method owners and decide, as their call, whether
      any criterion should be re-specified. All the evidence is in
      [`THRESHOLDS.md`](THRESHOLDS.md); none of it has been acted on:
  - `normalize_total` absolute tolerance is below one float32 ULP at the data's magnitude.
  - `umap.*.trustworthiness` scores one implementation, not agreement, and the CPU reference
    fails it on pbmc3k.
  - `umap.cross_embedding_knn_overlap` demands more agreement than reseeding the CPU gives.
  - `calculate_niche` UTAG thresholds (ARI/NMI 0.85) sit above the reference's own
    leidenalg-vs-igraph agreement of 0.5041, and `cluster_count_difference <= 1` is not
    achievable across Leiden implementations.
- [x] Probed the failures for an actual rapids-singlecell defect (jobs `38957332`, `38957359`,
      `38957415`) and found none. On the `calculate_niche` `neighborhood` kNN graph
      **rapids-singlecell is the accurate side**: it matches exact float64 ground truth to 4.9e-06
      while `sc.pp.neighbors` returns strictly worse neighbours on 743 of 4668 rows, because the
      profile holds only 755 distinct rows across 4668 cells. Leiden modularity, PCA subspace, and
      `normalize_total` were all cleared too. See the table in [`THRESHOLDS.md`](THRESHOLDS.md).
      This corrects an earlier note that read the divergence as an rsc-side defect.

## Upstream reports to file

None of these are fixed here. Each belongs to a different project.

- [ ] **scanpy — real bug, file first.** `sc.pp.neighbors` removes the cell itself from the k-NN
      result by position, assuming it sorts first. With exact duplicates it does not, so a genuine
      nearest neighbour is dropped and self keeps a slot: the row silently ends up with
      `n_neighbors - 2` real neighbours and distances shifted one place outward. Confirmed by exact
      containment on `squidpy.datasets.imc` (job `38957612`): 1261 rows retain self, 743 are inexact
      against float64 ground truth, all 743 are inside the 1261, and all 3407 rows that do not retain
      self are exact. Max excess distance +3.01. Reproducer, CPU-only and download-free:
      [`squidpy/scanpy_neighbors_duplicate_bug.py`](squidpy/scanpy_neighbors_duplicate_bug.py) —
      240 cells, 168 of them left with 13 real neighbours instead of 14. Versions: scanpy 1.12.3,
      anndata 0.13.2, pynndescent 0.6.0, scikit-learn 1.9.0.
- [ ] **rapids-singlecell — low-severity convention mismatch.** `rsc.pp.neighbors` materializes the
      self-loop as an explicit stored zero in `obsp["distances"]`, so it holds `n_neighbors` entries
      per row where scanpy holds `n_neighbors - 1` (700/700 rows on `pbmc68k_reduced`, nnz 10500 vs
      9800, ratio exactly 15/14). Neighbour sets are identical once the diagonal is removed and
      `obsp["connectivities"]` is unaffected (both 15900 nnz, no diagonal), so Leiden and UMAP are
      fine. Worth noting that this convention is what makes rsc immune to the scanpy bug above. Ask
      whether the extra diagonal entry is intended before treating it as a defect.
- [ ] **squidpy — method-design question, not a bug.** `calculate_niche(flavor="neighborhood")`
      clusters a feature space holding 755 distinct rows across 4668 cells, so 89.67% of cells have an
      exact distance tie at the k-th neighbour and the partition is not well determined. It also calls
      `sc.tl.leiden` without passing `flavor`, so results depend on which Leiden backend is installed —
      and those backends disagree substantially on this data (ARI 0.5041 between leidenalg and igraph
      on the UTAG space). Any CPU/GPU agreement threshold here measures tie-breaking and backend
      choice, not correctness.

## Cluster portability

- [ ] Resolve compatibility of the locked `rapids-singlecell-cu12==0.16.1` wheel with the
      cluster's V100 nodes. Jobs `38938670` and `38938725` reached CUDA but RSC preprocessing
      kernels failed with `named symbol not found`; A100 (`38936727`, `38938766`) and H100
      (`38956547`) nodes run the kernels. Prefer a wheel with the required architecture coverage
      so the checked-in job can remain a generic GPU request. The smoke check now detects this in
      seconds instead of part-way through the suite, but does not fix it.
- [ ] Report `gpusrv29` to the HPC team if it remains unhealthy. Job `38938173` saw the
      allocated A100 in `nvidia-smi`, but CUDA returned `cudaErrorNoDevice`; `nvidia-smi -q`
      reported an unknown MIG state and requested a GPU reset.

## Scientific validation

- [ ] Review every tolerance with the method owners and record a `basis` for each. All of them
      currently carry an empty basis, including the four diagnosed above.
- [ ] Decide, with the method owners, whether stochastic criteria should be expressed against a
      measured baseline rather than a hard-coded number. Evidence that a hard-coded number cannot
      travel: the reseeded-CPU embedding-overlap baseline moved from 0.556 to 0.586 between
      machines. This is a specification change and is not made unilaterally.
- [ ] Run stochastic comparisons over multiple fixed seeds and report distributions or
      confidence intervals, not only a single realization.
- [ ] Add at least one larger real dataset after the PBMC3k biological pipeline is stable,
      then report label-transfer accuracy, marker recovery, clustering agreement, and
      runtime/memory alongside low-level numerical metrics.
- [ ] Pin or checksum every externally downloaded dataset and record dataset provenance in
      the generated summary.
- [ ] Have an independent maintainer review the CPU reference choice and matched parameters
      for each method, especially Harmony and downstream biological metrics.

## Automation and reviewer evidence

- [ ] Configure a self-hosted GPU runner with the labels used by
      `.github/workflows/gpu-equivalence.yml`, then exercise the manual workflow and verify that
      reports and per-script logs upload even when thresholds fail.
- [ ] Add scheduled or release-triggered GPU runs once runner capacity and data access are
      reliable. Keep pull-request execution manual or scoped until the resource cost is understood.
- [ ] Link a stable generated `report/summary.md`, metric CSV, environment lock, GPU
      inventory, and raw logs from the reviewer response. State explicitly that automatic
      enforcement requires GPU CI.
- [ ] Decide which failures block releases: infrastructure and missing records should always
      block; provisional scientific thresholds should become blocking only after they are
      reviewed and baselined. The smoke check's exit 90 is the first infrastructure gate.

## Cluster notes

The Lustre checkout is a plain copy with no `.git`, so changes must be copied up before
submitting:

```bash
scp <changed files> hpc-submit01.scidom.de:/lustre/groups/ml01/workspace/$USER/rapids-singlecell-reproducibility/<path>/
```

Submit a full run, or a subset while iterating:

```bash
sbatch cluster/run_comparisons.sbatch
sbatch --time=01:00:00 --export=ALL,EQUIVALENCE_SCRIPTS='squidpy/calculate_niche.py' cluster/run_comparisons.sbatch
```

Results land in `benchmarks/comparison/cluster-results/<job-id>` in the shared checkout.
