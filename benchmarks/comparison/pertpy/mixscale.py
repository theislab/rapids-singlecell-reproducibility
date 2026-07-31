from __future__ import annotations

import pertpy as pt
import rapids_singlecell as rsc
from _report import lower, upper, write_report
from _shared import max_abs, pearson, screen_adata

if not hasattr(pt.tl, "Mixscale"):
    raise RuntimeError(
        "The installed pertpy version does not expose pt.tl.Mixscale; install a release containing Mixscale before generating equivalence evidence."
    )

cpu, gpu = screen_adata(), screen_adata()
cpu.layers["X_pert"] = cpu.X.copy()
gpu.layers["X_pert"] = gpu.X.copy()
pt.tl.Mixscale().mixscale(cpu, pert_key="gene_target", control="NT", layer="X_pert", test_method="t-test")
rsc.ptg.Mixscale().mixscale(gpu, pert_key="gene_target", control="NT", layer="X_pert", test_method="t-test")
cpu_score = cpu.obs["mixscale_score"].to_numpy(dtype=float)
gpu_score = gpu.obs["mixscale_score"].to_numpy(dtype=float)
metrics = [
    upper("mixscale.score_max_abs_error", max_abs(cpu_score, gpu_score), 1e-5),
    lower("mixscale.score_correlation", pearson(cpu_score, gpu_score), 0.9999),
]
write_report("mixscale", "seeded synthetic perturbation screen", "deterministic", metrics)
