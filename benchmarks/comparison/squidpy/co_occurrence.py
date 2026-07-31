from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import squidpy as sq

from _report import lower_bound, upper_bound, write_report
from _shared import (
    IMC_CLUSTER_KEY,
    load_imc,
    max_abs_error,
    mean_abs_error,
    pearson_correlation,
)


adata = load_imc()

# Use one explicit interval array for both implementations. This tests the statistic,
# rather than tiny differences in their independent min/max distance heuristics.
spatial = np.asarray(adata.obsm["spatial"], dtype=np.float32)
sample = spatial[: min(1000, adata.n_obs)]
distances = np.linalg.norm(sample[:, None, :] - sample[None, :, :], axis=2)
positive_distances = distances[distances > 0]
interval = np.linspace(
    np.quantile(positive_distances, 0.01),
    np.quantile(positive_distances, 0.50),
    50,
    dtype=np.float32,
)

reference_occ, reference_interval = sq.gr.co_occurrence(
    adata, cluster_key=IMC_CLUSTER_KEY, interval=interval, copy=True
)
candidate_occ, candidate_interval = rsc.gr.co_occurrence(
    adata, cluster_key=IMC_CLUSTER_KEY, interval=interval, copy=True
)

np.testing.assert_array_equal(reference_occ.shape, candidate_occ.shape)
metrics = [
    upper_bound(
        "interval.max_abs_error",
        max_abs_error(reference_interval, candidate_interval),
        1e-6,
    ),
    upper_bound(
        "occurrence.max_abs_error",
        max_abs_error(reference_occ, candidate_occ),
        1e-5,
    ),
    upper_bound(
        "occurrence.mean_abs_error",
        mean_abs_error(reference_occ, candidate_occ),
        1e-6,
    ),
    lower_bound(
        "occurrence.pearson_correlation",
        pearson_correlation(reference_occ, candidate_occ),
        0.99999,
    ),
]

write_report(
    method="co_occurrence",
    dataset="squidpy.datasets.imc",
    tier="deterministic",
    metrics=metrics,
)
