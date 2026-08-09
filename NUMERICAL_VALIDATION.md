# Numerical validation: the stated standard against the measured result

For assessment alongside the manuscript's **Methods — Numerical validation**.

The Methods make a testable claim with no free parameters:

> Deterministic numerical operations (normalization, HVG selection, PCA) produce outputs within
> floating-point tolerance (`numpy.allclose` with default parameters).

Applied literally — elementwise, `rtol=1e-5`, `atol=1e-8` — that criterion **holds for normalization
and PCA, and fails for scaling, regression, Pearson residuals, normalized dispersions and two
activity-inference outputs**. Nine of 47 such comparisons fail. The sentence needs qualifying before
publication.

## How it was measured

`numpy.allclose` accepts a pair of arrays when every element satisfies
`|a - b| <= atol + rtol * |b|`. The suite reports that as one scale-free number, `allclose_excess`:
the worst elementwise difference divided by its own envelope. At or below `1` the arrays are
`allclose`; above it they are not, and the value says by how much.

Two properties of the criterion matter for what follows. It is **relative**, whereas the suite
previously used hand-picked absolute tolerances that appear nowhere in the manuscript. And it is
**elementwise**, so a single element decides it.

Each failing criterion also records the reference magnitude at its worst element, which says which of
the two terms decided the verdict rather than leaving it to be inferred.

**These figures are one dated measurement, not a live claim.** They come from a single complete
20-group run: NVIDIA H100 80GB, rapids-singlecell 0.16.1, scanpy 1.12.3, squidpy 1.8.3, pertpy 1.1.1,
decoupler 2.2.0, 202 gating metrics of which 47 are `allclose` comparisons. This document is the only
place the assessment is written down; the current state of any run is in its own generated
`report/summary.md`, which lists every criterion with the measured diagnosis behind each failure. If
the two disagree, the report is right and this document is out of date.

## Result: 9 of 47 fail, for two different reasons

| Comparison                         | Excess | \|b\| at worst element | Decided by | Pearson, same arrays |
| ---------------------------------- | -----: | ---------------------: | ---------- | -------------------: |
| `regress_out`                      |   52.5 |                1.4e-06 | `atol`     |           1.00000000 |
| `normalize_pearson_residuals`      |   30.6 |                3.5e-04 | `atol`     |           1.00000000 |
| `hvg.seurat.dispersions_norm`      |   5.90 |                4.3e-04 | `atol`     |                    — |
| `perturbation_signature`           |   3.27 |                1.6e-04 | `atol`     |           1.00000000 |
| `scale`                            |   1.65 |                7.2e-04 | `atol`     |           1.00000000 |
| `ulm.adjusted_pvalue`              |   21.2 |                  0.998 | **`rtol`** |           1.00000000 |
| `mlm.score`                        |   15.4 |                 0.0109 | **`rtol`** |           1.00000000 |
| `hvg.cell_ranger.dispersions_norm` |   13.1 |                1.4e-03 | **`rtol`** |                    — |
| `mlm.adjusted_pvalue`              |   1.82 |                  0.979 | **`rtol`** |           1.00000000 |

**Five are an artefact of `atol`.** `atol=1e-8` is an absolute floor calibrated for float64. Where a
quantity passes through zero — residuals, z-scores, normalized dispersions, signed differences — that
floor rather than the relative term decides the criterion, and float32 cannot resolve a difference
below it there. No implementation can satisfy the default parameters on these quantities.

**Four are genuine relative disagreements.** At |b| ≈ 0.98–1.0 the absolute floor is irrelevant, so
`ulm.adjusted_pvalue` and `mlm.adjusted_pvalue` are real differences of roughly 2e-4 and 2e-5 relative
— on the order of 10² to 10³ float32 ULPs, larger than rounding alone. These are not artefacts of the
criterion. They are small in absolute terms and do not move any downstream result, but they exceed the
tolerance the manuscript cites.

## Per operation, against the sentence as written

| Operation named or implied in Methods              | Holds? | Worst excess |
| -------------------------------------------------- | ------ | -----------: |
| Normalization — `normalize_total`, `log1p`, `sqrt` | yes    |        0.015 |
| PCA — `explained_variance_ratio`                   | yes    |      5.2e-10 |
| HVG selection — means, variances                   | yes    |      1.6e-08 |
| HVG selection — `dispersions_norm`                 | **no** |         13.1 |
| Scaling — `scale`                                  | **no** |         1.65 |
| Regression — `regress_out`                         | **no** |         52.5 |
| Pearson residuals                                  | **no** |         30.6 |
| Activity inference — scores, adjusted p-values     | **no** |         21.2 |

The Methods name _normalization, HVG selection, PCA_. Normalization and PCA hold comfortably. HVG
selection holds for means and variances but not for normalized dispersions, so the claim is true for
part of a named operation and false for another part of it. The abstract's broader statement that
preprocessing functions were "confirmed numerically equivalent within floating-point tolerance"
covers scaling, which does not hold under the stated standard.

## This is not a rapids-singlecell defect

The suite measures a Pearson correlation on the same arrays that feed each `allclose` check, and for
every failure with a paired correlation it is `1.00000000`, against thresholds of 0.999 to 0.99999.
The implementations agree about the data; a small number of elements decide the criterion.

The same effect, seen from the other end of the scale, is what previously made `normalize_total` look
like a disagreement: an absolute tolerance of 1e-5 applied to values reaching 1751.05, where one
float32 ULP is already 1.220703125e-4. Under the manuscript's own relative criterion it passes with
143x headroom.

## The values are not stable across GPU architectures

The same suite, same image, same seeds, run on an A100 (sm80, driver 12.9) and an H100 (sm90, driver
13.0):

| Comparison                    | A100  | H100  |  Ratio |
| ----------------------------- | ----- | ----- | -----: |
| `normalize_total`             | 0.015 | 0.015 |  1.00x |
| `normalize_pearson_residuals` | 30.63 | 30.63 |  1.00x |
| `ulm.adjusted_pvalue`         | 21.18 | 21.18 |  1.00x |
| `regress_out`                 | 67.44 | 52.53 |  0.78x |
| `hvg.seurat.dispersions_norm` | 3.464 | 5.899 |  1.70x |
| `ulm.score`                   | 0.626 | 0.250 |  0.40x |
| `mlm.score`                   | 3.895 | 15.44 | **4x** |

Some deterministic-tier comparisons reproduce to the last digit across architectures; others move by
up to a factor of four. This is consistent with the manuscript's own caveat about accumulation order
and non-deterministic atomics, but it means a single reported tolerance for these operations is
hardware-specific. It is **not** a claim about same-hardware reproducibility, which was not tested.

## What is and is not being claimed

- **Not a defect report.** Correlation `1.00000000` on every failing array with a paired correlation
  says the GPU and CPU implementations agree.
- **No threshold was widened to reach these numbers.** The failures are recorded as failures and carry
  their measured diagnosis in the generated report. See [`OPEN.md`](OPEN.md) for the suite's other
  limits.
- **Half the failures are the criterion's fault, half are not.** Conflating them would be the easy
  mistake in either direction.
- **Replacement wording is deliberately not proposed here.** What the Methods should say is an
  authors' decision, and it depends on which standard reviewers are asked to accept.

The choice is between three: state the defaults as written and acknowledge the operations that do not
meet them; state a relative-only criterion appropriate to float32; or make a per-operation statement
saying which operations are exact, which agree within a stated envelope, and which are equivalent in
distribution rather than elementwise.

**The relative-only option does not work, at least for `regress_out`.** With the raw outputs now
stored, this is checkable rather than assumed. Dropping `atol` makes that comparison _worse_, not
better — worst excess rises from 167 to 6480 — because residuals that sit near zero have a large
relative error precisely where the absolute floor was protecting them. The single worst element is
`atol`-dominated, but only 31.8% of the 1,828 violating elements are, so "an artefact of `atol`" is
not the whole story for this operation either.

What the distribution does show is that the two implementations agree almost everywhere: 1,828 of
5,400,000 elements (0.034%) fall outside the envelope, the median excess is 0.021, and the 99.99th
percentile is 2.18. The criterion is decided by a very thin tail.

This was measured from the stored arrays on a laptop, with no GPU and no rerun. The same check has
not yet been done for the other eight failures.

## Reproducing this

```bash
docker run --rm --gpus all -v "$PWD/out:/out" rsc-equivalence
```

Every value above comes from that run's own machine-readable output — `out/equivalence.json` and
`out/report/metrics.csv` — and no criterion was altered to produce it.
