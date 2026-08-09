"""How two arrays are compared. One copy, applied at evaluation time.

These used to live in four near-identical `_shared.py` modules and run on the GPU node,
which meant the definitions could drift apart and a new comparison needed a new run. They
are pure functions of a stored reference/candidate pair, so they belong here.

`reference` is always the CPU side. Several of these are asymmetric — `allclose_excess`
builds its envelope from the reference — so the argument order matters.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

RTOL = 1e-5
ATOL = 1e-8


def dense(value) -> np.ndarray:
    if sparse.issparse(value):
        return value.toarray()
    return np.asarray(value)


def _pair(candidate, reference) -> tuple[np.ndarray, np.ndarray]:
    return (
        dense(candidate).astype(float, copy=False),
        dense(reference).astype(float, copy=False),
    )


def allclose_excess(candidate, reference, *, rtol: float = RTOL, atol: float = ATOL) -> float:
    """Worst elementwise violation of `numpy.allclose`, as a fraction of its own envelope.

    `numpy.allclose` accepts when every element satisfies |a - b| <= atol + rtol * |b|.
    Dividing by that envelope gives one scale-free number: `<= 1` means the arrays are
    `allclose`, and the value says how far past it the worst element sits. NaN propagates,
    because NaN is not close to anything and dropping it would make this weaker than the
    criterion it claims to be.
    """
    a, b = _pair(candidate, reference)
    return float(np.max(np.abs(a - b) / (atol + rtol * np.abs(b))))


def allclose_worst_magnitude(candidate, reference, *, rtol: float = RTOL, atol: float = ATOL) -> float:
    """Reference magnitude at the element that decides `allclose_excess`.

    This is what says which of the criterion's two terms actually bound it. Below
    `atol / rtol` the absolute floor dominates, and at float32 precision no implementation
    can satisfy it there; above, the disagreement is genuinely relative.
    """
    a, b = _pair(candidate, reference)
    ratio = (np.abs(a - b) / (atol + rtol * np.abs(b))).ravel()
    return float(np.abs(b.ravel()[int(np.nanargmax(ratio))]))


def max_abs_error(candidate, reference) -> float:
    a, b = _pair(candidate, reference)
    return float(np.nanmax(np.abs(a - b)))


def mean_abs_error(candidate, reference) -> float:
    a, b = _pair(candidate, reference)
    return float(np.nanmean(np.abs(a - b)))


def max_rel_error(candidate, reference) -> float:
    """Largest absolute difference relative to the reference's own scale."""
    a, b = _pair(candidate, reference)
    scale = float(np.nanmax(np.abs(b)))
    if scale == 0.0:
        return 0.0
    return float(np.nanmax(np.abs(a - b)) / scale)


def float32_ulp_at_max(candidate, reference) -> float:
    """One float32 ULP at the largest magnitude in the reference."""
    _, b = _pair(candidate, reference)
    return float(np.spacing(np.float32(float(np.nanmax(np.abs(b))))))


def pearson_correlation(candidate, reference) -> float:
    a, b = _pair(candidate, reference)
    a, b = a.ravel(), b.ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 2 or np.std(a[mask]) == 0 or np.std(b[mask]) == 0:
        return float(a[mask].shape == b[mask].shape and np.allclose(a[mask], b[mask]))
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def standard_deviation_max_abs_error(candidate, reference) -> float:
    a, b = _pair(candidate, reference)
    return float(np.nanmax(np.abs(a.std(axis=0) - b.std(axis=0))))


# Metric suffix -> how to compute it from a captured pair. A criterion named
# `<point>.<suffix>` is computed here whenever `<point>` was captured; anything absent
# falls back to a scalar the measurement script recorded, so a committed snapshot with no
# stored arrays still evaluates.
COMPARISONS = {
    "allclose_excess": allclose_excess,
    "pairwise_allclose_excess": allclose_excess,
    "score_allclose_excess": allclose_excess,
    "allclose_worst_magnitude": allclose_worst_magnitude,
    "max_abs_error": max_abs_error,
    "mean_abs_error": mean_abs_error,
    "max_rel_error": max_rel_error,
    "float32_ulp_at_max": float32_ulp_at_max,
    "pearson_correlation": pearson_correlation,
    "pairwise_correlation": pearson_correlation,
    "score_correlation": pearson_correlation,
    "standard_deviation_max_abs_error": standard_deviation_max_abs_error,
}
