from __future__ import annotations

import numpy as np
import pertpy as pt
import rapids_singlecell as rsc
from _report import measure, write_report
from _shared import allclose_excess, grouped_adata

adata = grouped_adata()
metrics = []
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
    excess = allclose_excess(candidate.to_numpy(), reference.to_numpy())
    metrics.append(measure(f"{name}.pairwise_allclose_excess", excess))
    metrics.append(
        measure(
            f"{name}.pairwise_correlation",
            np.corrcoef(reference.to_numpy().ravel(), candidate.to_numpy().ravel())[0, 1],
        )
    )

write_report("distance", "seeded grouped Gaussian data", "deterministic", metrics)
