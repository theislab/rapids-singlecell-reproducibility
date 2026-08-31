from __future__ import annotations

import numpy as np
import pertpy as pt
import rapids_singlecell as rsc
from _report import capture, write_report
from _shared import grouped_adata

METHOD = "distance"
# The reference group `onesided_distances` measures from. Any group works; g0 is the first.
SELECTED_GROUP = "g0"

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


def aligned(reference, candidate):
    """Both sides as 1-D arrays in the reference's own group order.

    `onesided_distances` returns a Series on one side and may return a single-column
    DataFrame on the other, and neither promises an order. Reindexing rather than
    comparing positionally means a mismatch shows up as a missing group, not as a
    numerical disagreement.
    """
    reference = reference.squeeze()
    candidate = candidate.squeeze().reindex(reference.index)
    assert not candidate.isna().any(), f"groups missing from the GPU result: {candidate[candidate.isna()].index}"
    return reference.to_numpy(), candidate.to_numpy()


for name in metric_names:
    reference = pt.tl.Distance(name, obsm_key="X_pca").pairwise(adata, groupby="group", show_progressbar=False)
    candidate = rsc.ptg.Distance(name, obsm_key="X_pca").pairwise(adata, groupby="group", multi_gpu=False)
    np.testing.assert_array_equal(reference.index, candidate.index)
    np.testing.assert_array_equal(reference.columns, candidate.columns)
    capture(METHOD, f"{name}.pairwise", reference=reference.to_numpy(), candidate=candidate.to_numpy())

    # `onesided_distances` is a separate code path, not a slice of the pairwise matrix: it
    # builds only the cross-group and diagonal pairs it needs.
    reference_one = pt.tl.Distance(name, obsm_key="X_pca").onesided_distances(
        adata, groupby="group", selected_group=SELECTED_GROUP, show_progressbar=False
    )
    candidate_one = rsc.ptg.Distance(name, obsm_key="X_pca").onesided_distances(
        adata, groupby="group", selected_group=SELECTED_GROUP, multi_gpu=False
    )
    reference_one, candidate_one = aligned(reference_one, candidate_one)
    capture(METHOD, f"{name}.onesided", reference=reference_one, candidate=candidate_one)

write_report(METHOD, "seeded grouped Gaussian data", "deterministic", [], shape=adata.shape)
