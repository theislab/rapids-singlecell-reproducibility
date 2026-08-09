from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import squidpy as sq
from _report import measure, write_report
from _shared import (
    IMC_CLUSTER_KEY,
    allclose_excess,
    load_imc,
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
    measure("interval.allclose_excess", allclose_excess(candidate_interval, reference_interval)),
    measure("occurrence.allclose_excess", allclose_excess(candidate_occ, reference_occ)),
    measure("occurrence.mean_abs_error", mean_abs_error(reference_occ, candidate_occ)),
    measure("occurrence.pearson_correlation", pearson_correlation(reference_occ, candidate_occ)),
]

write_report(
    method="co_occurrence",
    dataset="squidpy.datasets.imc",
    tier="deterministic",
    metrics=metrics,
)
