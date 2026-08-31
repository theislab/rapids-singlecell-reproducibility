from __future__ import annotations

import numpy as np
import pertpy as pt
import rapids_singlecell as rsc
from _report import capture, write_report
from _shared import guide_adata

METHOD = "guide_assignment"

metrics = []
for method in ("assign_by_threshold", "assign_to_max_guide"):
    cpu, gpu = guide_adata(), guide_adata()
    kwargs = {"assignment_threshold": 5}
    getattr(pt.pp.GuideAssignment(), method)(cpu, **kwargs)
    getattr(rsc.ptg.GuideAssignment(), method)(gpu, **kwargs)
    if method == "assign_by_threshold":
        cpu_assignment = cpu.layers["assigned_guides"]
        gpu_assignment = rsc.get.X_to_CPU(gpu.layers["assigned_guides"])
        if hasattr(cpu_assignment, "toarray"):
            cpu_assignment = cpu_assignment.toarray()
        if hasattr(gpu_assignment, "toarray"):
            gpu_assignment = gpu_assignment.toarray()
        capture(METHOD, method, reference=np.asarray(cpu_assignment), candidate=np.asarray(gpu_assignment))
    else:
        capture(METHOD, method, reference=cpu.obs["assigned_guide"], candidate=gpu.obs["assigned_guide"])

cpu, gpu = guide_adata(), guide_adata()
pt.pp.GuideAssignment().assign_mixture_model(cpu)
rsc.ptg.GuideAssignment().assign_mixture_model(gpu)
capture(METHOD, "assign_mixture_model", reference=cpu.obs["assigned_guide"], candidate=gpu.obs["assigned_guide"])

write_report(METHOD, "seeded Poisson guide-count mixture", "near-deterministic", metrics, shape=cpu.shape)
