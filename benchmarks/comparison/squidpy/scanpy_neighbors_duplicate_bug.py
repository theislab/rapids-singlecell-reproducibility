"""Minimal reproducer: sc.pp.neighbors loses a neighbour on rows with exact duplicates.

CPU-only, no GPU, no downloads, ~240 cells, deterministic. Written to be pasted into an
upstream scanpy issue.

`sc.pp.neighbors(n_neighbors=k)` requests k nearest neighbours (the cell itself included)
and stores k-1 per row in `obsp["distances"]`, removing self **by position** on the
assumption that self sorts first. When a cell has exact duplicates that assumption breaks:
a genuine nearest neighbour is dropped instead, self keeps a slot, and the row is silently
left with k-2 real neighbours and a distance vector shifted one place outward.

Observed with scanpy 1.12.3 / anndata 0.13.2 / pynndescent 0.6.0 / scikit-learn 1.9.0.

Real-world impact: on `squidpy.datasets.imc` the neighborhood-profile feature space used by
`sq.gr.calculate_niche(flavor="neighborhood")` holds 755 distinct rows across 4668 cells.
1261 rows retain self, and 743 of them return neighbour distances up to 3.01 larger than
exact float64 k-NN. Any downstream Leiden or UMAP consumes that graph. Duplicate rows are
routine in compositional, categorical, binned, or low-dimensional feature spaces.

Run:
    python benchmarks/comparison/squidpy/scanpy_neighbors_duplicate_bug.py
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

import numpy as np
import scanpy as sc
from anndata import AnnData

N_NEIGHBORS = 15
N_PROTOTYPES = 12
COPIES = 20


def main() -> None:
    for package in ("scanpy", "anndata", "pynndescent", "scikit-learn"):
        try:
            print(f"{package}=={version(package)}")
        except PackageNotFoundError:
            print(f"{package}==not installed")

    # Every cell has exactly COPIES-1 duplicates, so self rarely sorts first.
    rng = np.random.default_rng(0)
    prototypes = rng.normal(size=(N_PROTOTYPES, 4)).astype(np.float32)
    features = np.repeat(prototypes, COPIES, axis=0)
    n = features.shape[0]

    adata = AnnData(features.copy())
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X", random_state=0)
    distances = adata.obsp["distances"].tocsr()

    stored = np.diff(distances.indptr)
    retains_self = np.array(
        [i in set(distances.indices[distances.indptr[i] : distances.indptr[i + 1]].tolist()) for i in range(n)]
    )
    real = stored - retains_self.astype(int)
    expected = N_NEIGHBORS - 1

    print(f"\n{n} cells from {N_PROTOTYPES} distinct rows, each repeated {COPIES} times")
    print(f"n_neighbors={N_NEIGHBORS}, so every row should hold {expected} real neighbours")
    print(f"  stored entries per row : {sorted(set(stored.tolist()))}")
    print(f"  rows retaining self    : {int(retains_self.sum())}/{n}")
    print(f"  real neighbours per row: {dict(zip(*np.unique(real, return_counts=True), strict=True))}")
    short = int((real < expected).sum())
    print(f"\n  rows with fewer than {expected} real neighbours: {short}/{n}")
    assert short > 0, "expected the duplicate handling to drop a neighbour; it did not"
    print(f"  -> {short} rows silently carry {expected - 1} real neighbours instead of {expected}")


if __name__ == "__main__":
    main()
