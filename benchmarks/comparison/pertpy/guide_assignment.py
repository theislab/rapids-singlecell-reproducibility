from __future__ import annotations

import numpy as np
import pertpy as pt
import rapids_singlecell as rsc
from _report import measure, write_report
from _shared import guide_adata

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
        agreement = np.mean(np.asarray(cpu_assignment) == np.asarray(gpu_assignment))
    else:
        agreement = np.mean(cpu.obs["assigned_guide"].to_numpy() == gpu.obs["assigned_guide"].to_numpy())
    metrics.append(measure(f"{method}.assignment_agreement", agreement))

cpu, gpu = guide_adata(), guide_adata()
pt.pp.GuideAssignment().assign_mixture_model(cpu)
rsc.ptg.GuideAssignment().assign_mixture_model(gpu)
metrics.append(
    measure(
        "assign_mixture_model.assignment_agreement",
        np.mean(cpu.obs["assigned_guide"].to_numpy() == gpu.obs["assigned_guide"].to_numpy()),
    )
)

write_report("guide_assignment", "seeded Poisson guide-count mixture", "near-deterministic", metrics)
