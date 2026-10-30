from __future__ import annotations

import numpy as np
import rapids_singlecell as rsc
import squidpy as sq

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report
from benchmarks.comparison.squidpy._shared import load_imc

METHOD = "spatial_autocorr"

adata = load_imc()
# Float64 avoids treating expected float32 round-off as an algorithmic difference.
adata.X = adata.X.astype(np.float64)

for mode, statistic in (("moran", "I"), ("geary", "C")):
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

    np.testing.assert_array_equal(reference.index, candidate.index)
    capture(METHOD, f"{mode}.{statistic}", reference=reference[statistic], candidate=candidate[statistic])

write_report(
    method=METHOD,
    dataset="squidpy.datasets.imc",
    metrics=[],
    shape=adata.shape,
    reference_package="squidpy",
)
