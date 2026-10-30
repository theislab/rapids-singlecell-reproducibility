from __future__ import annotations

import anndata as ad
import numpy as np
import pertpy as pt
import rapids_singlecell as rsc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report

METHOD = "guide_assignment"


def guide_counts():
    rng = np.random.default_rng(42)
    x = np.zeros((200, 8), dtype=np.float32)
    for guide in range(8):
        x[:, guide] = np.concatenate([rng.poisson(2, 140), rng.poisson(50, 60)])
        rng.shuffle(x[:, guide])
    return ad.AnnData(x)


cpu, gpu = guide_counts(), guide_counts()
pt.pp.GuideAssignment().assign_by_threshold(cpu, assignment_threshold=5)
rsc.ptg.GuideAssignment().assign_by_threshold(gpu, assignment_threshold=5)
capture(METHOD, "assign_by_threshold", reference=cpu.layers["assigned_guides"], candidate=gpu.layers["assigned_guides"])

cpu, gpu = guide_counts(), guide_counts()
pt.pp.GuideAssignment().assign_to_max_guide(cpu, assignment_threshold=5)
rsc.ptg.GuideAssignment().assign_to_max_guide(gpu, assignment_threshold=5)
capture(METHOD, "assign_to_max_guide", reference=cpu.obs["assigned_guide"], candidate=gpu.obs["assigned_guide"])

cpu, gpu = guide_counts(), guide_counts()
pt.pp.GuideAssignment().assign_mixture_model(cpu)
rsc.ptg.GuideAssignment().assign_mixture_model(gpu)
capture(METHOD, "assign_mixture_model", reference=cpu.obs["assigned_guide"], candidate=gpu.obs["assigned_guide"])

write_report(
    METHOD,
    "seeded Poisson guide-count mixture",
    [],
    shape=cpu.shape,
    reference_package="pertpy",
)
