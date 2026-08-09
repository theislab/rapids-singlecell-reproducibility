"""Check that every copy of `allclose_excess` agrees with `numpy.allclose` itself.

The comparison packages each carry their own helper module, so this criterion exists in
several places. A copy that drifts would silently change the suite's verdict, so all of
them are checked against the function they are meant to reproduce.

    python benchmarks/comparison/test_allclose_excess.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
MODULES = [
    HERE / "scanpy_core" / "_shared.py",
    HERE / "scanpy_extended" / "_shared.py",
    HERE / "pertpy" / "_shared.py",
    HERE / "squidpy" / "_shared.py",
    HERE / "decoupler" / "decoupler.py",
]


def load(path: Path, name: str = "allclose_excess"):
    """Import one helper module without running its package's comparison scripts."""
    if path.name == "decoupler.py":
        # Importing this module runs the whole decoupler comparison, so read the
        # function out of the source instead of executing the script.
        namespace: dict = {"np": np}
        source = path.read_text()
        start = source.index(f"def {name}(")
        end = source.index("\n\n\n", start)
        exec(compile(source[start:end], str(path), "exec"), namespace)
        return namespace[name]
    spec = importlib.util.spec_from_file_location(f"helper_{path.parent.name}", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return getattr(module, name)


def main() -> int:
    rng = np.random.default_rng(0)
    cases = [
        # Identical inputs, at magnitudes where an absolute tolerance would disagree.
        (np.zeros(64), np.zeros(64)),
        (np.full(64, 1751.05), np.full(64, 1751.05)),
        # One float32 ULP apart at a large magnitude: allclose, absolute error 1.2e-4.
        (np.float32(1751.05) + np.float32(np.spacing(np.float32(1751.05))), np.float32(1751.05)),
        # Differences straddling the envelope, and small values where atol dominates.
        (rng.normal(size=256), rng.normal(size=256)),
        (rng.normal(size=256) * 1e-9, rng.normal(size=256) * 1e-9),
        (np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0 + 1e-3])),
    ]
    for path in MODULES:
        excess = load(path)
        for left, right in cases:
            left_array = np.atleast_1d(np.asarray(left, dtype=float))
            right_array = np.atleast_1d(np.asarray(right, dtype=float))
            value = excess(left_array, right_array)
            expected = np.allclose(left_array, right_array)
            assert (value <= 1.0) == expected, f"{path.parent.name}: excess={value} but np.allclose={expected}"
            assert value >= 0.0, f"{path.parent.name}: negative excess {value}"
        print(f"ok  {path.relative_to(HERE)}")

    # The diagnosis has to name the term that actually decided the verdict.
    for path in (
        HERE / "scanpy_core" / "_shared.py",
        HERE / "pertpy" / "_shared.py",
        HERE / "decoupler" / "decoupler.py",
    ):
        diagnose = load(path, "allclose_diagnosis")
        near_zero = diagnose(np.array([1.0, 1e-9]), np.array([1.0, 3e-9]))
        assert "absolute floor" in near_zero, f"{path.parent.name}: atol case not identified: {near_zero}"
        well_scaled = diagnose(np.array([1.0, 500.0]), np.array([1.0, 400.0]))
        assert "relative disagreement" in well_scaled, f"{path.parent.name}: rtol case not identified: {well_scaled}"
        print(f"ok  {path.relative_to(HERE)} (diagnosis)")

    print(f"\nAll {len(MODULES)} copies agree with numpy.allclose on {len(cases)} cases.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
