from __future__ import annotations

import pertpy as pt
import rapids_singlecell as rsc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.pertpy._shared import screen_adata
from benchmarks.comparison.report import write_report

METHOD = "mixscape"

cpu, gpu = screen_adata(), screen_adata()
pt.tl.Mixscape().perturbation_signature(cpu, pert_key="gene_target", control="NT", n_neighbors=5)
rsc.ptg.Mixscape().perturbation_signature(
    gpu, pert_key="gene_target", control="NT", n_neighbors=5, knn_algorithm="brute"
)
capture(METHOD, "perturbation_signature", reference=cpu.layers["X_pert"], candidate=gpu.layers["X_pert"])

pt.tl.Mixscape().mixscape(cpu, pert_key="gene_target", control="NT", test_method="t-test")
rsc.ptg.Mixscape().mixscape(gpu, pert_key="gene_target", control="NT", test_method="t-test")
capture(
    METHOD,
    "mixscape.global_class",
    reference=cpu.obs["mixscape_class_global"],
    candidate=gpu.obs["mixscape_class_global"],
)
capture(METHOD, "mixscape.p_ko", reference=cpu.obs["mixscape_class_p_ko"], candidate=gpu.obs["mixscape_class_p_ko"])

pt.tl.Mixscape().lda(cpu, pert_key="gene_target", control="NT", test_method="t-test")
rsc.ptg.Mixscape().lda(gpu, pert_key="gene_target", control="NT", test_method="t-test")
capture(METHOD, "mixscape.lda", reference=cpu.uns["mixscape_lda"], candidate=gpu.uns["mixscape_lda"])

write_report(
    METHOD,
    "seeded synthetic perturbation screen",
    [],
    shape=cpu.shape,
    reference_package="pertpy",
)
