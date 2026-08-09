from __future__ import annotations

from itertools import combinations

import cupy as cp
import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
import squidpy as sq
from _report import capture, write_report
from _shared import dataframe_values

METHOD = "ligrec"

adata = sc.datasets.paul15()
clusters = adata.obs["paul15_clusters"].cat.categories[:4].tolist()
adata = adata[adata.obs["paul15_clusters"].isin(clusters)].copy()
adata.obs["paul15_clusters"] = adata.obs["paul15_clusters"].cat.remove_unused_categories()
sc.pp.normalize_total(adata)
adata.raw = adata.copy()

# Fix both the interaction set and cluster pairs so the comparison never depends on
# an external OmniPath query and remains small enough for routine reproduction.
genes = adata.var_names[:8].tolist()
interactions = list(combinations(genes, 2))
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
np.testing.assert_array_equal(reference["pvalues"].columns, candidate["pvalues"].columns)

reference_means = dataframe_values(reference["means"])
candidate_means = dataframe_values(candidate["means"])
reference_pvalues = dataframe_values(reference["pvalues"])
candidate_pvalues = dataframe_values(candidate["pvalues"])

# p-values legitimately carry NaN where a test could not be run, so the NaN pattern is
# compared as its own criterion and the numeric criteria run over the finite entries.
capture(METHOD, "means", reference=reference_means, candidate=candidate_means)
capture(METHOD, "pvalues", reference=reference_pvalues, candidate=candidate_pvalues)

write_report(
    method=METHOD,
    dataset="scanpy.datasets.paul15",
    tier="stochastic",
    metrics=[],
)
