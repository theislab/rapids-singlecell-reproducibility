from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import squidpy as sq

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report
from benchmarks.comparison.squidpy._shared import IMC_CLUSTER_KEY, load_imc

METHOD = "co_occurrence"

adata = load_imc()

# Use one explicit interval array for both implementations. This tests the statistic,
# rather than tiny differences in their independent min/max distance heuristics.
spatial = np.asarray(adata.obsm["spatial"], dtype=np.float32)
sample = spatial[:1000]
distances = np.linalg.norm(sample[:, None, :] - sample[None, :, :], axis=2)
positive_distances = distances[distances > 0]
interval = np.linspace(
    np.quantile(positive_distances, 0.01),
    np.quantile(positive_distances, 0.50),
    50,
    dtype=np.float32,
)

reference, _ = sq.gr.co_occurrence(adata, cluster_key=IMC_CLUSTER_KEY, interval=interval, copy=True)
candidate, _ = rsc.gr.co_occurrence(adata, cluster_key=IMC_CLUSTER_KEY, interval=interval, copy=True)

capture(METHOD, "occurrence", reference=reference, candidate=candidate)

write_report(
    method=METHOD,
    dataset="squidpy.datasets.imc",
    metrics=[],
    shape=adata.shape,
    reference_package="squidpy",
)
