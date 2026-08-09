from __future__ import annotations

import inspect

import numpy as np
import rapids_singlecell as rsc
import scanpy as sc
from _report import measure, write_report
from _shared import graph_overlap, pearson

metrics = []
adata = sc.datasets.pbmc68k_reduced()
adata.obs["batch"] = np.where(np.arange(adata.n_obs) % 2, "batch_a", "batch_b")
reference, candidate = adata.copy(), adata.copy()
sc.external.pp.bbknn(reference, batch_key="batch", neighbors_within_batch=3, n_pcs=50)
rsc.pp.bbknn(candidate, batch_key="batch", neighbors_within_batch=3, n_pcs=50, algorithm="brute")
metrics.append(
    measure(
        "bbknn.connectivity_jaccard", graph_overlap(reference.obsp["connectivities"], candidate.obsp["connectivities"])
    )
)

counts = sc.datasets.pbmc3k()
reference, candidate = counts.copy(), counts.copy()
rsc.get.anndata_to_GPU(candidate)
cpu_seed_arg = {"rng": 0} if "rng" in inspect.signature(sc.pp.scrublet).parameters else {"random_state": 0}
sc.pp.scrublet(reference, n_prin_comps=30, use_approx_neighbors=True, **cpu_seed_arg)
rsc.pp.scrublet(candidate, random_state=0, n_prin_comps=30, use_approx_neighbors=True, verbose=False)
metrics.extend(
    [
        measure("scrublet.score_correlation", pearson(reference.obs["doublet_score"], candidate.obs["doublet_score"])),
        measure(
            "scrublet.call_agreement",
            np.mean(reference.obs["predicted_doublet"].to_numpy() == candidate.obs["predicted_doublet"].to_numpy()),
        ),
    ]
)

reference_sim = sc.pp.scrublet_simulate_doublets(counts, sim_doublet_ratio=1.0, random_seed=0)
counts_gpu = counts.copy()
rsc.get.anndata_to_GPU(counts_gpu)
candidate_sim = rsc.pp.scrublet_simulate_doublets(counts_gpu, sim_doublet_ratio=1.0, random_seed=0)
ref_totals = np.asarray(reference_sim.X.sum(axis=1)).ravel()
candidate_sim_x = rsc.get.X_to_CPU(candidate_sim.X)
gpu_totals = np.asarray(candidate_sim_x.sum(axis=1)).ravel()
metrics.extend(
    [
        measure(
            "scrublet_simulate_doublets.mean_library_size_relative_error",
            abs(ref_totals.mean() - gpu_totals.mean()) / ref_totals.mean(),
        ),
        measure(
            "scrublet_simulate_doublets.std_library_size_relative_error",
            abs(ref_totals.std() - gpu_totals.std()) / ref_totals.std(),
        ),
    ]
)

write_report("bbknn_scrublet", "pbmc68k_reduced + pbmc3k", "stochastic", metrics)
