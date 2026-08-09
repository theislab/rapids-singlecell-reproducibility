from __future__ import annotations

import numpy as np
import pertpy as pt
import rapids_singlecell as rsc
from _report import measure, write_report
from _shared import allclose_excess, pearson, screen_adata

cpu, gpu = screen_adata(), screen_adata()
pt.tl.Mixscape().perturbation_signature(cpu, pert_key="gene_target", control="NT", n_neighbors=5)
rsc.ptg.Mixscape().perturbation_signature(
    gpu, pert_key="gene_target", control="NT", n_neighbors=5, knn_algorithm="brute"
)
cpu_signature = cpu.layers["X_pert"].toarray() if hasattr(cpu.layers["X_pert"], "toarray") else cpu.layers["X_pert"]
gpu_signature = rsc.get.X_to_CPU(gpu.layers["X_pert"])
if hasattr(gpu_signature, "toarray"):
    gpu_signature = gpu_signature.toarray()
metrics = [
    measure("perturbation_signature.allclose_excess", allclose_excess(gpu_signature, cpu_signature)),
    measure("perturbation_signature.pearson_correlation", pearson(cpu_signature, gpu_signature)),
]

pt.tl.Mixscape().mixscape(cpu, pert_key="gene_target", control="NT", test_method="t-test")
rsc.ptg.Mixscape().mixscape(gpu, pert_key="gene_target", control="NT", test_method="t-test")
metrics.extend(
    [
        measure(
            "mixscape.global_class_agreement",
            np.mean(cpu.obs["mixscape_class_global"].to_numpy() == gpu.obs["mixscape_class_global"].to_numpy()),
        ),
        measure("mixscape.p_ko_correlation", pearson(cpu.obs["mixscape_class_p_ko"], gpu.obs["mixscape_class_p_ko"])),
    ]
)

pt.tl.Mixscape().lda(cpu, pert_key="gene_target", control="NT", test_method="t-test")
rsc.ptg.Mixscape().lda(gpu, pert_key="gene_target", control="NT", test_method="t-test")
metrics.append(measure("mixscape.lda_abs_correlation", abs(pearson(cpu.uns["mixscape_lda"], gpu.uns["mixscape_lda"]))))

write_report("mixscape", "seeded synthetic perturbation screen", "near-deterministic", metrics)
