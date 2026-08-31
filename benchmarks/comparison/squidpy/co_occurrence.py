from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import squidpy as sq
from _report import capture, write_report
from _shared import IMC_CLUSTER_KEY, load_imc

METHOD = "co_occurrence"

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
capture(METHOD, "interval", reference=reference_interval, candidate=candidate_interval)
capture(METHOD, "occurrence", reference=reference_occ, candidate=candidate_occ)

write_report(
    method=METHOD,
    dataset="squidpy.datasets.imc",
    tier="deterministic",
    metrics=[],
    shape=adata.shape,
)
