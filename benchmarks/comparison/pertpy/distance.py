from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pertpy as pt
import rapids_singlecell as rsc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report

METHOD = "distance"
rng = np.random.default_rng(42)
x = np.vstack([rng.normal(loc=i * 0.5, size=(30, 12)) for i in range(4)])
adata = ad.AnnData(x, obs=pd.DataFrame({"group": pd.Categorical([f"g{i}" for i in range(4) for _ in range(30)])}))
adata.obsm["X_pca"] = x.copy()
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

    reference_one = pt.tl.Distance(name, obsm_key="X_pca").onesided_distances(
        adata, groupby="group", selected_group="g0", show_progressbar=False
    )
    candidate_one = rsc.ptg.Distance(name, obsm_key="X_pca").onesided_distances(
        adata, groupby="group", selected_group="g0", multi_gpu=False
    )
    reference_one, candidate_one = reference_one.squeeze(), candidate_one.squeeze()
    np.testing.assert_array_equal(reference_one.index, candidate_one.index)
    capture(METHOD, f"{name}.onesided", reference=reference_one.to_numpy(), candidate=candidate_one.to_numpy())

write_report(METHOD, "seeded grouped Gaussian data", [], shape=adata.shape, reference_package="pertpy")
