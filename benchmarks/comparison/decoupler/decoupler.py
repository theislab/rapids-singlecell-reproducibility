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

sys.path.insert(0, str(Path(__file__).parents[1]))

import decoupler as dc
import numpy as np
import rapids_singlecell as rsc
from arrays import capture

METHOD = "decoupler_methods"


def package_version(package):
    candidates = (package, "rapids-singlecell-cu12") if package == "rapids-singlecell" else (package,)
    for candidate in candidates:
        try:
            return version(candidate)
        except PackageNotFoundError:
            pass
    return "unknown"


def compare_frames(name, reference, candidate):
    """Store the two frames; evaluate.py compares them."""
    np.testing.assert_array_equal(reference.index, candidate.index)
    np.testing.assert_array_equal(reference.columns, candidate.columns)
    capture(METHOD, name, reference=reference.to_numpy(dtype=float), candidate=candidate.to_numpy(dtype=float))
    return []


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
    "method": METHOD,
    "reference_package": "decoupler",
    "dataset": "decoupler.ds.toy",
    "tier": "deterministic",
    "versions": {"rapids-singlecell": package_version("rapids-singlecell"), "decoupler": package_version("decoupler")},
    "metrics": metrics,
}
output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "decoupler_methods.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
# Measuring is not evaluating: the verdict is formed by evaluate.py from this record.
