"""Self-check for the parts of the evaluation path that decide a verdict.

Not a test suite for the science — the comparisons are checked by the run itself. This
covers only the logic that could turn a bad measurement into a green criterion, which no
GPU run would reveal because it needs a degenerate input to trigger.

    python benchmarks/comparison/test_evaluation.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import scipy.sparse as sp
from comparisons import (
    abs_allclose_excess,
    allclose_excess,
    graph_exact_agreement,
    graph_jaccard,
    pearson_correlation,
    relative_l2_max,
)
from evaluate import decide, enriched


def test_all_nan_candidate_does_not_correlate() -> None:
    """An absent measurement must not read as perfect agreement."""
    reference = np.arange(100.0)
    assert math.isnan(pearson_correlation(np.full(100, np.nan), reference))
    assert math.isnan(pearson_correlation(np.full(100, np.nan), np.full(100, np.nan)))
    assert not decide({"metric": "x", "observed": float("nan"), "comparison": ">=", "tolerance": 0.999})
    # A constant pair is still answerable: undefined correlation, but they do agree.
    assert pearson_correlation(np.ones(10), np.ones(10)) == 1.0
    assert pearson_correlation(np.zeros(10), np.ones(10)) == 0.0
    assert abs(pearson_correlation(reference, reference) - 1.0) < 1e-12


def test_strict_comparison_is_not_silently_relaxed() -> None:
    """`>` used to fall through to `>=`, which would gate at a threshold nobody chose."""
    at_threshold = {"metric": "harmony", "observed": 0.95, "tolerance": 0.95}
    assert decide({**at_threshold, "comparison": ">="})
    assert not decide({**at_threshold, "comparison": ">"})
    assert decide({**at_threshold, "observed": 0.951, "comparison": ">"})
    try:
        decide({**at_threshold, "comparison": "~="})
    except SystemExit:
        pass
    else:
        raise AssertionError("an unknown comparison must not be evaluated as >=")


def test_enriched_record_keeps_measurements_and_drops_criteria() -> None:
    """A promoted record carries what was measured, never the threshold it was judged by."""
    record = {
        "method": "m",
        "versions": {},
        "metrics": [
            {"metric": "a.allclose_excess", "observed": 0.5, "tolerance": 1.0, "gating": True, "passed": True},
        ],
        "derived_metrics": ["a.allclose_excess"],
        "passed": True,
    }
    assert enriched(record) == {
        "method": "m",
        "versions": {},
        "metrics": [{"metric": "a.allclose_excess", "observed": 0.5}],
        "passed": True,
    }


def test_abs_allclose_ignores_component_sign_but_not_real_error() -> None:
    """A principal component's sign is arbitrary; a difference in its magnitude is not."""
    reference = np.random.default_rng(0).normal(size=(200, 10))
    flipped = reference * np.where(np.arange(10) % 2 == 0, -1.0, 1.0)
    # The signed envelope reads a sign flip as total disagreement, which is why the gate on
    # PCA scores and loadings uses magnitudes.
    assert allclose_excess(flipped, reference) > 1.0
    assert abs_allclose_excess(flipped, reference) == 0.0
    perturbed = reference.copy()
    perturbed[0, 0] += 1.0
    assert abs_allclose_excess(perturbed, reference) > 1.0


def test_relative_l2_flags_a_scaled_component_and_a_degenerate_reference() -> None:
    reference = np.random.default_rng(1).normal(size=(200, 5))
    assert relative_l2_max(reference, reference) == 0.0
    scaled = reference.copy()
    scaled[:, 3] *= 1.5
    assert abs(relative_l2_max(scaled, reference) - 0.5) < 1e-9
    # A zero-norm reference component must not read as agreement: any difference against it
    # is infinitely relative, and returning 0 there would hide the whole component.
    degenerate = reference.copy()
    degenerate[:, 0] = 0.0
    assert math.isinf(relative_l2_max(reference, degenerate))
    assert relative_l2_max(degenerate, degenerate) == 0.0


def test_graph_exact_agreement_is_all_or_nothing_per_row() -> None:
    """`graph_jaccard` gives partial credit; this counts only rows that match completely."""
    dense = 1.0 - np.eye(4)
    same = sp.csr_matrix(dense)
    assert graph_exact_agreement(same, same) == 1.0
    changed = dense.copy()
    changed[0, 1] = 0.0
    differing = sp.csr_matrix(changed)
    assert graph_exact_agreement(same, differing) == 0.75
    assert graph_jaccard(same, differing) > graph_exact_agreement(same, differing)
    # Self-loops must not decide it: the two implementations disagree on storing the diagonal.
    with_diagonal = same.tolil()
    with_diagonal[0, 0] = 1.0
    assert graph_exact_agreement(with_diagonal.tocsr(), same) == 1.0


if __name__ == "__main__":
    for name, case in sorted(globals().items()):
        if name.startswith("test_"):
            case()
            print(f"ok {name}")
