from __future__ import annotations

import pertpy as pt
import rapids_singlecell as rsc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.pertpy._shared import screen_adata
from benchmarks.comparison.report import write_report

METHOD = "mixscale"

cpu, gpu = screen_adata(), screen_adata()
cpu.layers["X_pert"] = cpu.X.copy()
gpu.layers["X_pert"] = gpu.X.copy()
pt.tl.Mixscale().mixscale(cpu, pert_key="gene_target", control="NT", layer="X_pert", test_method="t-test")
rsc.ptg.Mixscale().mixscale(gpu, pert_key="gene_target", control="NT", layer="X_pert", test_method="t-test")
capture(METHOD, "mixscale.score", reference=cpu.obs["mixscale_score"], candidate=gpu.obs["mixscale_score"])
write_report(METHOD, "seeded synthetic perturbation screen", [], shape=cpu.shape, reference_package="pertpy")
