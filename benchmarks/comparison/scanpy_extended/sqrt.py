from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
from _report import upper, write_report

reference = sc.datasets.pbmc3k()
candidate = reference.copy()
rsc.get.anndata_to_GPU(candidate)
sc.pp.sqrt(reference)
rsc.pp.sqrt(candidate)
candidate_x = rsc.get.X_to_CPU(candidate.X)
if hasattr(reference.X, "toarray"):
    reference_x = reference.X.toarray()
    candidate_x = candidate_x.toarray()
else:
    reference_x = np.asarray(reference.X)
    candidate_x = np.asarray(candidate_x)
error = np.max(np.abs(reference_x - candidate_x))
write_report("sqrt", "scanpy.datasets.pbmc3k", "deterministic", [upper("X.max_abs_error", error, 1e-6)])
