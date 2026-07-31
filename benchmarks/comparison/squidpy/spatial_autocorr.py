from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import squidpy as sq
from _report import lower_bound, upper_bound, write_report
from _shared import load_imc, max_abs_error, pearson_correlation

adata = load_imc()
metrics = []

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

    metrics.extend(
        [
            upper_bound(
                f"{mode}.{statistic}.max_abs_error",
                max_abs_error(ref_values, candidate_values),
                1e-6,
            ),
            lower_bound(
                f"{mode}.{statistic}.pearson_correlation",
                pearson_correlation(ref_values, candidate_values),
                0.999999,
            ),
        ]
    )

write_report(
    method="spatial_autocorr",
    dataset="squidpy.datasets.imc",
    tier="deterministic",
    metrics=metrics,
)
