from __future__ import annotations

import pertpy as pt
import rapids_singlecell as rsc
from _report import capture, write_report
from _shared import screen_adata

METHOD = "mixscape"

cpu, gpu = screen_adata(), screen_adata()
pt.tl.Mixscape().perturbation_signature(cpu, pert_key="gene_target", control="NT", n_neighbors=5)
rsc.ptg.Mixscape().perturbation_signature(
    gpu, pert_key="gene_target", control="NT", n_neighbors=5, knn_algorithm="brute"
)
cpu_signature = cpu.layers["X_pert"].toarray() if hasattr(cpu.layers["X_pert"], "toarray") else cpu.layers["X_pert"]
gpu_signature = rsc.get.X_to_CPU(gpu.layers["X_pert"])
if hasattr(gpu_signature, "toarray"):
    gpu_signature = gpu_signature.toarray()
capture(METHOD, "perturbation_signature", reference=cpu_signature, candidate=gpu_signature)

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
# The LDA axis has no fixed sign, so the criterion is on the absolute correlation.
capture(METHOD, "mixscape.lda", reference=cpu.uns["mixscape_lda"], candidate=gpu.uns["mixscape_lda"])

write_report(METHOD, "seeded synthetic perturbation screen", "near-deterministic", [])
