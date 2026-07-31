from __future__ import annotations

from itertools import combinations

import cupy as cp
import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
import squidpy as sq

from _report import lower_bound, upper_bound, write_report
from _shared import (
    dataframe_values,
    max_abs_error,
    mean_abs_error,
    pearson_correlation,
)


adata = sc.datasets.paul15()
sc.pp.normalize_total(adata)
adata.raw = adata.copy()

# Fix both the interaction set and cluster pairs so the comparison never depends on
# an external OmniPath query and remains small enough for routine reproduction.
genes = adata.var_names[:8].tolist()
interactions = list(combinations(genes, 2))
clusters = adata.obs["paul15_clusters"].cat.categories[:4].tolist()
n_perms = 500

reference = sq.gr.ligrec(
    adata,
    "paul15_clusters",
    interactions=interactions,
    clusters=clusters,
    threshold=0.01,
    corr_method=None,
    use_raw=True,
    copy=True,
    n_perms=n_perms,
    seed=0,
    n_jobs=1,
    numba_parallel=False,
    show_progress_bar=False,
)

# rapids-singlecell currently seeds permutations through CuPy's global RNG.
cp.random.seed(0)
candidate = rsc.gr.ligrec(
    adata,
    "paul15_clusters",
    interactions=interactions,
    clusters=clusters,
    threshold=0.01,
    corr_method=None,
    use_raw=True,
    copy=True,
    n_perms=n_perms,
)

np.testing.assert_array_equal(reference["means"].index, candidate["means"].index)
np.testing.assert_array_equal(reference["means"].columns, candidate["means"].columns)
np.testing.assert_array_equal(reference["pvalues"].index, candidate["pvalues"].index)
np.testing.assert_array_equal(
    reference["pvalues"].columns, candidate["pvalues"].columns
)

reference_means = dataframe_values(reference["means"])
candidate_means = dataframe_values(candidate["means"])
reference_pvalues = dataframe_values(reference["pvalues"])
candidate_pvalues = dataframe_values(candidate["pvalues"])

valid = np.isfinite(reference_pvalues) & np.isfinite(candidate_pvalues)
nan_agreement = np.mean(np.isnan(reference_pvalues) == np.isnan(candidate_pvalues))

metrics = [
    upper_bound(
        "means.max_abs_error",
        max_abs_error(reference_means, candidate_means),
        1e-5,
    ),
    lower_bound("pvalues.nan_mask_agreement", nan_agreement, 1.0),
    upper_bound(
        "pvalues.mean_abs_error",
        mean_abs_error(reference_pvalues[valid], candidate_pvalues[valid]),
        0.05,
    ),
    lower_bound(
        "pvalues.pearson_correlation",
        pearson_correlation(reference_pvalues[valid], candidate_pvalues[valid]),
        0.90,
    ),
]

write_report(
    method="ligrec",
    dataset="scanpy.datasets.paul15",
    tier="stochastic",
    metrics=metrics,
)
