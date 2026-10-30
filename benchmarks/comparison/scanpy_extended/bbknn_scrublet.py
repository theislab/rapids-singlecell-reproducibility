from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report
from benchmarks.comparison.scanpy_core._shared import gpu_copy

METHOD = "bbknn_scrublet"

adata = sc.datasets.pbmc68k_reduced()
adata.obs["batch"] = np.where(np.arange(adata.n_obs) % 2, "batch_a", "batch_b")
reference, candidate = adata.copy(), adata.copy()
# approx=False: exact (cKDTree) neighbor search on both sides, matching rsc's algorithm="brute".
sc.external.pp.bbknn(reference, batch_key="batch", neighbors_within_batch=3, n_pcs=50, approx=False)
rsc.pp.bbknn(candidate, batch_key="batch", neighbors_within_batch=3, n_pcs=50, algorithm="brute")
capture(METHOD, "bbknn", reference=reference.obsp["connectivities"], candidate=candidate.obsp["connectivities"])

counts = sc.datasets.pbmc3k()
reference, candidate = counts.copy(), gpu_copy(counts)
sc.pp.scrublet(reference, n_prin_comps=30, use_approx_neighbors=True, rng=0)
rsc.pp.scrublet(candidate, random_state=0, n_prin_comps=30, use_approx_neighbors=True, verbose=False)
capture(METHOD, "scrublet.score", reference=reference.obs["doublet_score"], candidate=candidate.obs["doublet_score"])
capture(
    METHOD, "scrublet.call", reference=reference.obs["predicted_doublet"], candidate=candidate.obs["predicted_doublet"]
)

reference_sim = sc.pp.scrublet_simulate_doublets(counts, sim_doublet_ratio=1.0, rng=0)
candidate_sim = rsc.pp.scrublet_simulate_doublets(gpu_copy(counts), sim_doublet_ratio=1.0, rng=0)
# Doublets are random cell pairs, so only the library-size distribution can match.
capture(
    METHOD,
    "scrublet_simulate_doublets.library_size",
    reference=np.asarray(reference_sim.X.sum(axis=1)).ravel(),
    candidate=candidate_sim.X.sum(axis=1).ravel(),
)

write_report(METHOD, "pbmc68k_reduced + pbmc3k", [])
