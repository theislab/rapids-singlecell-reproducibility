from __future__ import annotations

import json
import os
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

# Avoid this script's filename shadowing the installed ``decoupler`` package.
script_dir = str(Path(__file__).parent)
if script_dir in sys.path:
    sys.path.remove(script_dir)

import decoupler as dc
import numpy as np
import rapids_singlecell as rsc


def package_version(package):
    candidates = (package, "rapids-singlecell-cu12") if package == "rapids-singlecell" else (package,)
    for candidate in candidates:
        try:
            return version(candidate)
        except PackageNotFoundError:
            pass
    return "unknown"


def metric(name, observed, comparison, tolerance, passed):
    return {
        "metric": name,
        "observed": float(observed),
        "comparison": comparison,
        "tolerance": float(tolerance),
        "criterion": f"{comparison} {tolerance}",
        "passed": bool(passed),
    }


def allclose_excess(candidate, reference, *, rtol=1e-5, atol=1e-8):
    """Worst elementwise violation of `numpy.allclose` at its default parameters.

    The published validation standard for deterministic operations is `numpy.allclose`,
    which is relative: |a - b| <= atol + rtol * |b|. Dividing by that envelope gives one
    scale-free number where <= 1 means the two arrays are `allclose`.
    """
    return float(np.max(np.abs(candidate - reference) / (atol + rtol * np.abs(reference))))


def allclose_diagnosis(candidate, reference, *, rtol=1e-5, atol=1e-8):
    """Explain a `numpy.allclose` verdict by locating where its envelope is spent.

    The criterion has two terms, |a - b| <= atol + rtol * |b|, and which one binds depends
    entirely on the magnitude of the element that fails. For a quantity that crosses zero
    the absolute floor binds, and at float32 precision that floor is below the noise of any
    implementation. Reporting the reference magnitude at the worst element says which of
    the two terms the verdict actually rests on.
    """
    left_array = np.asarray(candidate, dtype=float)
    right_array = np.asarray(reference, dtype=float)
    envelope = atol + rtol * np.abs(right_array)
    ratio = (np.abs(left_array - right_array) / envelope).ravel()
    magnitude = float(np.abs(right_array.ravel()[int(np.nanargmax(ratio))]))
    if atol >= rtol * magnitude:
        return (
            f"The worst element sits at a reference magnitude of {magnitude:.3g}, where "
            f"numpy.allclose's absolute floor (atol={atol:g}) decides the criterion instead of its "
            f"relative term. This quantity crosses zero, and float32 cannot resolve a difference "
            f"below that floor there, so no float32 implementation can meet the default parameters "
            f"on it. The paired correlation on the same arrays says how the values agree overall."
        )
    return (
        f"The worst element sits at a reference magnitude of {magnitude:.3g}, where the relative "
        f"term (rtol={rtol:g}) decides the criterion, so this is a genuine relative disagreement."
    )


def compare_frames(name, reference, candidate, *, min_correlation=0.999):
    np.testing.assert_array_equal(reference.index, candidate.index)
    np.testing.assert_array_equal(reference.columns, candidate.columns)
    a, b = reference.to_numpy(dtype=float), candidate.to_numpy(dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    excess = allclose_excess(b[mask], a[mask])
    correlation = 1.0 if np.array_equal(a[mask], b[mask]) else np.corrcoef(a[mask], b[mask])[0, 1]
    return [
        metric(f"{name}.allclose_excess", excess, "<=", 1.0, excess <= 1.0)
        | {"diagnosis": allclose_diagnosis(b[mask], a[mask])},
        metric(f"{name}.pearson_correlation", correlation, ">=", min_correlation, correlation >= min_correlation),
    ]


adata, net = dc.ds.toy(nobs=80, nvar=40, bval=2, seed=42, verbose=False)
net = dc.pp.prune(features=adata.var_names, net=net, tmin=3)
metrics = []

methods = {
    "mlm": {},
    "ulm": {},
    "aucell": {"n_up": 10},
    "zscore": {"flavor": "RoKAI"},
    # times=0 compares the deterministic weighted aggregation itself. Permutation
    # calibration is separately covered by the shared seeded implementation tests.
    "waggr": {"fun": "wmean", "times": 0, "seed": 42},
}

for name, kwargs in methods.items():
    cpu, gpu = adata.copy(), adata.copy()
    getattr(dc.mt, name)(cpu, net, tmin=3, verbose=False, **kwargs)
    getattr(rsc.dcg, name)(gpu, net, tmin=3, verbose=False, **kwargs)
    metrics.extend(compare_frames(f"{name}.score", cpu.obsm[f"score_{name}"], gpu.obsm[f"score_{name}"]))
    padj_key = f"padj_{name}"
    if padj_key in cpu.obsm and padj_key in gpu.obsm:
        metrics.extend(compare_frames(f"{name}.adjusted_pvalue", cpu.obsm[padj_key], gpu.obsm[padj_key]))

report = {
    "method": "decoupler_methods",
    "reference_package": "decoupler",
    "dataset": "decoupler.ds.toy",
    "tier": "deterministic",
    "versions": {"rapids-singlecell": package_version("rapids-singlecell"), "decoupler": package_version("decoupler")},
    "passed": all(item["passed"] for item in metrics),
    "metrics": metrics,
}
output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "decoupler_methods.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
if not report["passed"]:
    raise AssertionError("One or more Decoupler equivalence thresholds failed")
