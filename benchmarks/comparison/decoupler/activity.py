from __future__ import annotations

import decoupler as dc
import numpy as np
import rapids_singlecell as rsc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report

METHOD = "decoupler_methods"


def compare_frames(name, reference, candidate):
    np.testing.assert_array_equal(reference.index, candidate.index)
    np.testing.assert_array_equal(reference.columns, candidate.columns)
    capture(METHOD, name, reference=reference.to_numpy(dtype=float), candidate=candidate.to_numpy(dtype=float))


adata, net = dc.ds.toy(nobs=80, nvar=40, bval=2, seed=42, verbose=False)
net = dc.pp.prune(features=adata.var_names, net=net, tmin=3)

methods = {
    "mlm": {},
    "ulm": {},
    "aucell": {"n_up": 10},
    "zscore": {"flavor": "RoKAI"},
    "waggr": {"fun": "wmean", "times": 0},
}

for name, kwargs in methods.items():
    cpu, gpu = adata.copy(), adata.copy()
    getattr(dc.mt, name)(cpu, net, tmin=3, verbose=False, **kwargs)
    getattr(rsc.dcg, name)(gpu, net, tmin=3, verbose=False, **kwargs)
    compare_frames(f"{name}.score", cpu.obsm[f"score_{name}"], gpu.obsm[f"score_{name}"])
    padj_key = f"padj_{name}"
    if padj_key in cpu.obsm:
        compare_frames(f"{name}.adjusted_pvalue", cpu.obsm[padj_key], gpu.obsm[padj_key])

write_report(
    METHOD,
    "decoupler.ds.toy",
    [],
    reference_package="decoupler",
    shape=adata.shape,
)
