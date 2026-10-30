from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def screen_adata():
    rng = np.random.default_rng(1)
    x = np.hstack(
        [
            np.vstack([np.clip(rng.normal(0, 1, (10, 10)), 0, None) for _ in range(3)]),
            np.vstack(
                [
                    np.clip(rng.normal(4, 0.7, (10, 10)), 0, None),
                    np.clip(rng.normal(4, 0.7, (10, 10)), 0, None),
                    np.clip(rng.normal(7, 0.9, (10, 10)), 0, None),
                ]
            ),
        ]
    ).astype(np.float32)
    obs = pd.DataFrame({"gene_target": ["NT"] * 10 + ["target_gene_a"] * 20}, index=[str(i) for i in range(30)])
    var = pd.DataFrame(index=[f"gene{i}" for i in range(20)])
    return ad.AnnData(sparse.csr_matrix(x), obs=obs, var=var)
