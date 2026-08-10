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

from comparisons import pearson_correlation
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


if __name__ == "__main__":
    for name, case in sorted(globals().items()):
        if name.startswith("test_"):
            case()
            print(f"ok {name}")
