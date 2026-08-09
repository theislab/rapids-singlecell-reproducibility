from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import squidpy as sq
from _report import capture, write_report
from _shared import load_imc

METHOD = "spatial_autocorr"

adata = load_imc()

for mode, statistic in (("moran", "I"), ("geary", "C")):
    # Float64 avoids treating expected float32 round-off as an algorithmic difference.
    adata.X = adata.X.astype(np.float64)
    reference = sq.gr.spatial_autocorr(
        adata,
        mode=mode,
        n_perms=None,
        corr_method=None,
        copy=True,
        n_jobs=1,
        show_progress_bar=False,
    )
    candidate = rsc.gr.spatial_autocorr(
        adata,
        mode=mode,
        n_perms=None,
        corr_method=None,
        dtype=np.float64,
        copy=True,
    )

    common_genes = reference.index.intersection(candidate.index)
    np.testing.assert_array_equal(reference.index, candidate.index)
    ref_values = reference.loc[common_genes, statistic].to_numpy()
    candidate_values = candidate.loc[common_genes, statistic].to_numpy()

    capture(METHOD, f"{mode}.{statistic}", reference=ref_values, candidate=candidate_values)

write_report(
    method=METHOD,
    dataset="squidpy.datasets.imc",
    tier="deterministic",
    metrics=[],
)
