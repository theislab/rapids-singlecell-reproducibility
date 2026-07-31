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


def compare_frames(name, reference, candidate, *, max_error=1e-4, min_correlation=0.999):
    np.testing.assert_array_equal(reference.index, candidate.index)
    np.testing.assert_array_equal(reference.columns, candidate.columns)
    a, b = reference.to_numpy(dtype=float), candidate.to_numpy(dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    error = np.max(np.abs(a[mask] - b[mask]))
    correlation = 1.0 if np.array_equal(a[mask], b[mask]) else np.corrcoef(a[mask], b[mask])[0, 1]
    return [
        metric(f"{name}.max_abs_error", error, "<=", max_error, error <= max_error),
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
    score_error = 1e-3 if name == "aucell" else 1e-4
    metrics.extend(
        compare_frames(f"{name}.score", cpu.obsm[f"score_{name}"], gpu.obsm[f"score_{name}"], max_error=score_error)
    )
    padj_key = f"padj_{name}"
    if padj_key in cpu.obsm and padj_key in gpu.obsm:
        metrics.extend(
            compare_frames(f"{name}.adjusted_pvalue", cpu.obsm[padj_key], gpu.obsm[padj_key], max_error=5e-4)
        )

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
