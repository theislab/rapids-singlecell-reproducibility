from __future__ import annotations

import pertpy as pt
import rapids_singlecell as rsc
from _report import capture, write_report
from _shared import screen_adata

METHOD = "mixscale"

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
capture(METHOD, "mixscale.score", reference=cpu_score, candidate=gpu_score)
write_report(METHOD, "seeded synthetic perturbation screen", "deterministic", [], shape=cpu.shape)
