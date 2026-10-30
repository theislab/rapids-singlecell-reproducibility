from __future__ import annotations

from itertools import combinations

import cupy as cp
import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
import squidpy as sq

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report

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
    cluster_key="paul15_clusters",
    interactions=interactions,
    clusters=clusters,
    threshold=0.01,
    corr_method=None,
    use_raw=True,
    copy=True,
    n_perms=n_perms,
    rng=0,
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

for key in ("means", "pvalues"):
    np.testing.assert_array_equal(reference[key].index, candidate[key].index)
    np.testing.assert_array_equal(reference[key].columns, candidate[key].columns)
    capture(METHOD, key, reference=reference[key].to_numpy(dtype=float), candidate=candidate[key].to_numpy(dtype=float))

write_report(
    method=METHOD,
    dataset="scanpy.datasets.paul15",
    metrics=[],
    shape=adata.shape,
    reference_package="squidpy",
)
