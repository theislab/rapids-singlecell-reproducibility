from __future__ import annotations

import inspect

import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
from _report import capture, write_report

METHOD = "bbknn_scrublet"

adata = sc.datasets.pbmc68k_reduced()
adata.obs["batch"] = np.where(np.arange(adata.n_obs) % 2, "batch_a", "batch_b")
reference, candidate = adata.copy(), adata.copy()
sc.external.pp.bbknn(reference, batch_key="batch", neighbors_within_batch=3, n_pcs=50)
rsc.pp.bbknn(candidate, batch_key="batch", neighbors_within_batch=3, n_pcs=50, algorithm="brute")
capture(METHOD, "bbknn", reference=reference.obsp["connectivities"], candidate=candidate.obsp["connectivities"])

counts = sc.datasets.pbmc3k()
reference, candidate = counts.copy(), counts.copy()
rsc.get.anndata_to_GPU(candidate)
cpu_seed_arg = {"rng": 0} if "rng" in inspect.signature(sc.pp.scrublet).parameters else {"random_state": 0}
sc.pp.scrublet(reference, n_prin_comps=30, use_approx_neighbors=True, **cpu_seed_arg)
rsc.pp.scrublet(candidate, random_state=0, n_prin_comps=30, use_approx_neighbors=True, verbose=False)
capture(METHOD, "scrublet.score", reference=reference.obs["doublet_score"], candidate=candidate.obs["doublet_score"])
capture(
    METHOD, "scrublet.call", reference=reference.obs["predicted_doublet"], candidate=candidate.obs["predicted_doublet"]
)

reference_sim = sc.pp.scrublet_simulate_doublets(counts, sim_doublet_ratio=1.0, random_seed=0)
counts_gpu = counts.copy()
rsc.get.anndata_to_GPU(counts_gpu)
candidate_sim = rsc.pp.scrublet_simulate_doublets(counts_gpu, sim_doublet_ratio=1.0, random_seed=0)
ref_totals = np.asarray(reference_sim.X.sum(axis=1)).ravel()
candidate_sim_x = rsc.get.X_to_CPU(candidate_sim.X)
gpu_totals = np.asarray(candidate_sim_x.sum(axis=1)).ravel()
# Simulated doublets are drawn from random cell pairs, so the two sides cannot be compared
# cell by cell; the library-size distribution is what has to match.
capture(METHOD, "scrublet_simulate_doublets.library_size", reference=ref_totals, candidate=gpu_totals)

write_report(METHOD, "pbmc68k_reduced + pbmc3k", "stochastic", [])
