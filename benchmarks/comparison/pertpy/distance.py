from __future__ import annotations

import numpy as np
import pertpy as pt
import rapids_singlecell as rsc
from _report import capture, write_report
from _shared import grouped_adata

METHOD = "distance"

adata = grouped_adata()
metric_names = (
    "edistance",
    "euclidean",
    "root_mean_squared_error",
    "mse",
    "mean_absolute_error",
    "pearson_distance",
    "cosine_distance",
    "r2_distance",
    "wasserstein",
)

for name in metric_names:
    reference = pt.tl.Distance(name, obsm_key="X_pca").pairwise(adata, groupby="group", show_progressbar=False)
    candidate = rsc.ptg.Distance(name, obsm_key="X_pca").pairwise(adata, groupby="group", multi_gpu=False)
    np.testing.assert_array_equal(reference.index, candidate.index)
    np.testing.assert_array_equal(reference.columns, candidate.columns)
    capture(METHOD, f"{name}.pairwise", reference=reference.to_numpy(), candidate=candidate.to_numpy())

write_report(METHOD, "seeded grouped Gaussian data", "deterministic", [])
