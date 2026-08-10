# Numerical validation: the stated standard against the measured result

For assessment alongside the manuscript's **Methods — Numerical validation**.

The Methods make a testable claim with no free parameters:

> Deterministic numerical operations (normalization, HVG selection, PCA) produce outputs within
> floating-point tolerance (`numpy.allclose` with default parameters).

Applied literally — elementwise, `rtol=1e-5`, `atol=1e-8` — that criterion **holds for normalization
and PCA, and fails for scaling, regression, Pearson residuals, normalized dispersions and two
activity-inference outputs**. Nine of 50 such comparisons fail. The sentence needs qualifying before
publication.

The suite gates on this standard only where the publication asserts it — the operations the Methods
and abstract name. Activity inference, perturbation signatures and the pertpy distances are measured
against it and reported, but carry no verdict, because a threshold there would be the suite's
invention rather than the paper's claim. Both kinds appear below: what follows is an assessment of
the sentence, not a pass list.

## How it was measured

`numpy.allclose` accepts a pair of arrays when every element satisfies
`|a - b| <= atol + rtol * |b|`. The suite reports that as one scale-free number, `allclose_excess`:
the worst elementwise difference divided by its own envelope. At or below `1` the arrays are
`allclose`; above it they are not, and the value says by how much.

Two properties of the criterion matter. It is **relative**, whereas the suite previously used
hand-picked absolute tolerances that appear nowhere in the manuscript. And it is **elementwise**, so
a single element decides it.

The raw CPU and GPU outputs are stored, and every quantity below is computed from them at evaluation
time rather than fixed while the run was happening. That is what makes the last two columns of the
table possible: the reference magnitude at the deciding element, which says which of the criterion's
two terms bound it, and the fraction of elements actually outside the envelope.

**Every figure below is committed and checkable.** The run is
[`snapshots/2026-08-09-derived`](benchmarks/comparison/snapshots/2026-08-09-derived), whose records
re-score to the same verdicts with no GPU:

```bash
python benchmarks/comparison/evaluate.py --results benchmarks/comparison/snapshots/2026-08-09-derived/results
```

**These figures are one dated measurement, not a live claim.** Complete 20-group run on an NVIDIA
A100-PCIE-40GB (MIG 3g.20gb, driver 12.9), rapids-singlecell 0.16.1, scanpy 1.12.3, squidpy 1.8.3,
pertpy 1.1.1, decoupler 2.2.0. 161 gating metrics; 50 comparisons are measured against `allclose`,
of which 30 gate and 20 are recorded as evidence on operations the paper does not name. The
current state of any run is in its own generated `report/summary.md`; if the two disagree, the
report is right and this document is out of date.

## Result: 9 of 50 fail

| Comparison                         | Gates? | Excess | \|b\| at deciding element | Elements outside envelope | Decided by | Pearson, same arrays |
| ---------------------------------- | ------ | -----: | ------------------------: | ------------------------: | ---------- | -------------------: |
| `regress_out`                      | yes    |  67.44 |                  1.22e-04 |                   0.0367% | `atol`     |           1.00000000 |
| `normalize_pearson_residuals`      | yes    |  30.63 |                  3.51e-04 |                  0.00150% | `atol`     |           1.00000000 |
| `ulm.adjusted_pvalue`              | no     |  21.18 |                     0.998 |                    0.250% | **`rtol`** |           1.00000000 |
| `hvg.cell_ranger.dispersions_norm` | yes    |  13.07 |                  1.44e-03 |                    0.751% | **`rtol`** |                    — |
| `mlm.score`                        | no     |  3.895 |                  1.09e-02 |                    0.500% | **`rtol`** |           1.00000000 |
| `hvg.seurat.dispersions_norm`      | yes    |  3.464 |                  3.22e-04 |                    0.102% | `atol`     |                    — |
| `perturbation_signature`           | no     |  2.686 |                  7.41e-03 |                    1.000% | **`rtol`** |           1.00000000 |
| `mlm.adjusted_pvalue`              | no     |  1.821 |                     0.979 |                    0.750% | **`rtol`** |           1.00000000 |
| `scale`                            | yes    |  1.787 |                  7.23e-04 |                 0.000130% | `atol`     |           1.00000000 |

The four marked *no* are the operations the publication states no tolerance for. They are measured
because the alternative is worse: those method groups are otherwise covered by a correlation floor,
and correlation reads `1.00000000` on the very arrays sitting 21x outside the `allclose` envelope. The
number is kept; only the verdict is dropped.

**Every failure is decided by a thin tail.** The largest violating fraction is 1% and most are far
below that; `scale` fails on roughly one element in 770,000. Where a paired correlation exists it is
`1.00000000`. The two implementations agree about the data — a handful of elements decide the
criterion.

**Four are an artefact of `atol`.** `atol=1e-8` is an absolute floor calibrated for float64. Where a
quantity passes through zero — residuals, z-scores, normalized dispersions — that floor rather than
the relative term decides the criterion, and float32 cannot resolve a difference below it there. No
implementation can satisfy the default parameters on those.

**Five are genuine relative disagreements.** At `|b| ≈ 0.98–1.0` the absolute floor is irrelevant, so
`ulm.adjusted_pvalue` and `mlm.adjusted_pvalue` are real differences of roughly 2e-4 and 2e-5
relative — on the order of 10² to 10³ float32 ULPs, larger than rounding alone. They are small in
absolute terms and move no downstream result, but they exceed the tolerance the manuscript cites.

**The atol/rtol attribution is itself not stable.** Which term decides depends on which element
happens to be worst, and that element moves between runs — `perturbation_signature` has been seen
decided both ways. The split above describes this run, not a fixed property of each operation.

## Per operation, against the sentence as written

| Operation named or implied in Methods              | Holds? | Worst excess |
| -------------------------------------------------- | ------ | -----------: |
| Normalization — `normalize_total`, `log1p`, `sqrt` | yes    |        0.015 |
| PCA — `explained_variance_ratio`                   | yes    |      5.2e-10 |
| HVG selection — means, variances                   | yes    |        0.354 |
| HVG selection — `dispersions_norm`                 | **no** |        13.07 |
| Scaling — `scale`                                  | **no** |        1.787 |
| Regression — `regress_out`                         | **no** |        67.44 |
| Pearson residuals                                  | **no** |        30.63 |
| Activity inference — scores, adjusted p-values     | **no** |        21.18 |

The Methods name _normalization, HVG selection, PCA_. Normalization and PCA hold comfortably. HVG
selection holds for means and variances but not for normalized dispersions, so the claim is true for
part of a named operation and false for another part of it. The abstract's broader statement that
preprocessing functions were "confirmed numerically equivalent within floating-point tolerance"
covers scaling, which does not hold under the stated standard.

The nearest miss among the comparisons inside the envelope is `wasserstein` at 0.986 — within 1.5% of
exceeding it.

## This is not a rapids-singlecell defect

Where a paired Pearson correlation is measured on the same arrays it is `1.00000000`, against
thresholds of 0.999 to 0.99999, and the violating fractions above are all at or below 1%.

The same effect, seen from the other end of the scale, is what previously made `normalize_total` look
like a disagreement: an absolute tolerance of 1e-5 applied to values reaching 1751.05, where one
float32 ULP is already 1.220703125e-4. Under the manuscript's own relative criterion it passes with
143x headroom.

## The values are architecture-specific

Same suite, same image, same seeds:

| Comparison                    | A100 (two nodes) |  H100 |
| ----------------------------- | ---------------: | ----: |
| `normalize_total`             |            0.015 | 0.015 |
| `normalize_pearson_residuals` |            30.63 | 30.63 |
| `ulm.adjusted_pvalue`         |            21.18 | 21.18 |
| `scale`                       |            1.787 | 1.648 |
| `regress_out`                 |            67.44 | 52.53 |
| `mlm.score`                   |            3.895 | 15.44 |

Two separate A100 nodes reproduce these to the digit; the H100 does not. So the differences are a
property of the **architecture**, not run-to-run noise, and a single reported tolerance for these
operations is hardware-specific. This is consistent with the manuscript's own caveat about
accumulation order and non-deterministic atomics.

## What is and is not being claimed

- **Not a defect report.** Perfect correlation on every failing array with a paired correlation says
  the GPU and CPU implementations agree.
- **No threshold was widened to reach these numbers.** The failures are recorded as failures and
  carry their measured diagnosis in the generated report. See [`OPEN.md`](OPEN.md) for the suite's
  other limits.
- **Roughly half the failures are the criterion's fault and half are not.** Conflating them would be
  the easy mistake in either direction.
- **Replacement wording is deliberately not proposed here.** What the Methods should say is an
  authors' decision, and it depends on which standard reviewers are asked to accept.

The choice is between three: state the defaults as written and acknowledge the operations that do not
meet them; state a relative-only criterion appropriate to float32; or make a per-operation statement
saying which operations are exact, which agree within a stated envelope, and which are equivalent in
distribution rather than elementwise.

**The relative-only option does not work, at least for `regress_out`.** Because the raw outputs are
stored, this is checkable rather than assumed. Dropping `atol` makes that comparison _worse_, not
better — worst excess rises from 167 to 6480 on the run where it was checked — because residuals
sitting near zero have a large relative error precisely where the absolute floor was protecting them.
The same check has not yet been run for the other eight.

## Reproducing this

```bash
docker run --rm --gpus all -v "$PWD/out:/out" rsc-equivalence
```

Every value above comes from that run's own output and is committed in the snapshot, so a reader can
re-derive the table without a GPU. Asking a _different_ question of the same run — a different
metric, the shape of the error distribution — needs the raw Zarr stores, which stay with the run
rather than in the repository.
