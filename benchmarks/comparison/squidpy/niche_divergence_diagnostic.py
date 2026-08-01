"""Diagnostic for the calculate_niche CPU/GPU divergence. Not part of the gated suite.

`calculate_niche.py` records that the `neighborhood` and `utag` flavors disagree between
Squidpy and rapids-singlecell while `cellcharter` agrees closely. This script locates *where*
that disagreement originates instead of widening the tolerances, by separating three possible
sources:

1. the features fed to clustering (neighborhood profile, UTAG product, PCA),
2. the kNN graph built from those features, and
3. the Leiden implementation that partitions the graph.

For the graph it also asks which side is actually *right*, using exact float64 kNN as the
arbiter — graph disagreement on its own says nothing about correctness.

It prints measurements and asserts nothing, so it never gates a run. Run it directly:

    python benchmarks/comparison/squidpy/niche_divergence_diagnostic.py

Findings from jobs 38956728, 38956767, 38957332, 38957359 and 38957415 are written up in
../THRESHOLDS.md. Summary: no rapids-singlecell defect. On this feature space rsc's kNN is
exact while scanpy's is not, because the profile holds 3913 duplicate rows out of 4668.
"""

from __future__ import annotations

import traceback

import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
import squidpy as sq
from _shared import IMC_CLUSTER_KEY, load_imc
from anndata import AnnData
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors

RESOLUTION = 0.5
N_NEIGHBORS = 15


def banner(text: str) -> None:
    print(f"\n{'=' * 78}\n{text}\n{'=' * 78}", flush=True)


def describe(name: str, labels) -> None:
    counts = pd.Series(np.asarray(labels).astype(str)).value_counts()
    print(
        f"{name:28} n_clusters={counts.size:3d}  "
        f"sizes: max={counts.iloc[0]:5d} median={int(counts.median()):5d} min={counts.iloc[-1]:5d}"
    )


def agreement(name: str, left, right) -> None:
    left = np.asarray(left).astype(str)
    right = np.asarray(right).astype(str)
    print(
        f"{name:44} ARI={adjusted_rand_score(left, right):.6f}  "
        f"NMI={normalized_mutual_info_score(left, right):.6f}  "
        f"n_clusters {pd.unique(left).size} vs {pd.unique(right).size}"
    )


def graph_jaccard(left, right) -> float:
    left = left.tocsr()
    right = right.tocsr()
    scores = []
    for index in range(left.shape[0]):
        a = set(left.indices[left.indptr[index] : left.indptr[index + 1]]) - {index}
        b = set(right.indices[right.indptr[index] : right.indptr[index + 1]]) - {index}
        scores.append(len(a & b) / max(1, len(a | b)))
    return float(np.mean(scores))


def scaled_neighborhood_profile(adata) -> np.ndarray:
    from squidpy.gr._niche import _calculate_neighborhood_profile

    matrix = adata.obsp["spatial_connectivities"].tocoo()
    profile = _calculate_neighborhood_profile(adata, IMC_CLUSTER_KEY, matrix, False)
    values = np.asarray(profile.to_numpy() if hasattr(profile, "to_numpy") else profile, dtype=np.float64)
    return ((values - values.mean(axis=0)) / (values.std(axis=0) + 1e-12)).astype(np.float32)


def utag_feature_matrix(adata) -> np.ndarray:
    adjacency = adata.obsp["spatial_connectivities"]
    row_sums = np.asarray(adjacency.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    normalized = sparse.diags(1.0 / row_sums) @ adjacency
    dense_x = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    return np.asarray(normalized @ dense_x).astype(np.float32)


def strip_diagonal(matrix):
    """Drop the diagonal, preserving explicitly stored zeros elsewhere.

    `nonzero()` hides stored zeros, and this feature space is full of genuine zero
    distances between duplicate cells, so the structure is walked directly.
    """
    matrix = matrix.tocsr()
    data, indices, indptr = [], [], [0]
    for row in range(matrix.shape[0]):
        span = slice(matrix.indptr[row], matrix.indptr[row + 1])
        for index, value in zip(matrix.indices[span], matrix.data[span], strict=True):
            if index != row:
                indices.append(index)
                data.append(value)
        indptr.append(len(indices))
    return sparse.csr_matrix((np.asarray(data), np.asarray(indices), np.asarray(indptr)), shape=matrix.shape)


def accuracy_against_exact(features, *, width: int = 13) -> None:
    """Which implementation actually returns the nearest neighbors?

    Graph disagreement alone says nothing about correctness. Exact float64 kNN is the
    arbiter: a correct implementation can never return a neighbor distance smaller than
    the true one, so entries materially larger than exact are genuine misses.
    """
    features = np.ascontiguousarray(features)
    exact = NearestNeighbors(n_neighbors=N_NEIGHBORS + 1, algorithm="brute", metric="euclidean")
    exact.fit(features.astype(np.float64))
    true_distances, _ = exact.kneighbors(features.astype(np.float64))
    reference = true_distances[:, 1 : width + 1]

    cpu = AnnData(features.copy())
    sc.pp.neighbors(cpu, n_neighbors=N_NEIGHBORS, use_rep="X", random_state=0)
    gpu = AnnData(features.copy())
    rsc.get.anndata_to_GPU(gpu)
    rsc.pp.neighbors(gpu, n_neighbors=N_NEIGHBORS, use_rep="X", algorithm="brute")
    rsc.get.anndata_to_CPU(gpu)

    print(f"\naccuracy vs exact float64 kNN, first {width} non-self neighbors per row:")
    for label, matrix in (("cpu", cpu.obsp["distances"]), ("gpu", gpu.obsp["distances"])):
        stripped = strip_diagonal(matrix).tocsr()
        candidate = np.full((stripped.shape[0], width), np.inf)
        for row in range(stripped.shape[0]):
            values = np.sort(stripped.data[stripped.indptr[row] : stripped.indptr[row + 1]])
            candidate[row, : min(width, values.size)] = values[: min(width, values.size)]
        delta = candidate - reference
        worse = delta > 1e-4 * np.maximum(1.0, np.abs(reference))
        print(
            f"  {label}: entries worse than exact={worse.sum():6d}  rows affected={worse.any(axis=1).sum():5d}  "
            f"max excess={delta.max():+.3e}"
        )

    # Convention difference worth knowing about: rsc materializes the self-loop as a
    # stored zero in `distances`, so it holds one more entry per row than scanpy. The
    # real neighbor sets match once the diagonal is removed, and `connectivities` is
    # unaffected.
    for key in ("distances", "connectivities"):
        cpu_matrix, gpu_matrix = cpu.obsp[key].tocsr(), gpu.obsp[key].tocsr()
        cpu_diagonal = sum(
            1
            for i in range(cpu_matrix.shape[0])
            if i in set(cpu_matrix.indices[cpu_matrix.indptr[i] : cpu_matrix.indptr[i + 1]].tolist())
        )
        gpu_diagonal = sum(
            1
            for i in range(gpu_matrix.shape[0])
            if i in set(gpu_matrix.indices[gpu_matrix.indptr[i] : gpu_matrix.indptr[i + 1]].tolist())
        )
        print(
            f"  obsp['{key}']: nnz cpu={cpu_matrix.nnz} gpu={gpu_matrix.nnz}; "
            f"rows with a stored diagonal cpu={cpu_diagonal} gpu={gpu_diagonal}"
        )


def compare_leiden_backends(features, tag: str) -> None:
    """Partition one identical graph with every Leiden backend.

    The neighbor graph is built once on the CPU and reused, so any label difference comes
    from community detection alone. Scanpy's two flavors bound how much the reference
    implementation disagrees with itself.
    """
    reference = AnnData(np.asarray(features, dtype=np.float32))
    sc.pp.neighbors(reference, n_neighbors=N_NEIGHBORS, use_rep="X", random_state=0)

    labels = {}
    for flavor in ("leidenalg", "igraph"):
        candidate = reference.copy()
        kwargs = {"flavor": flavor, "key_added": "labels", "resolution": RESOLUTION, "random_state": 0}
        if flavor == "igraph":
            kwargs["n_iterations"] = 2
        try:
            sc.tl.leiden(candidate, **kwargs)
            labels[f"cpu:{flavor}"] = candidate.obs["labels"].to_numpy()
        except Exception:  # noqa: BLE001 - a missing backend must not stop the comparison
            print(f"  cpu:{flavor} failed:\n{traceback.format_exc()}")

    gpu = reference.copy()
    try:
        rsc.get.anndata_to_GPU(gpu)
        rsc.tl.leiden(gpu, resolution=RESOLUTION, random_state=0, key_added="labels")
        labels["gpu:cugraph"] = np.asarray(gpu.obs["labels"])
    except Exception:  # noqa: BLE001 - a missing backend must not stop the comparison
        print(f"  gpu:cugraph failed:\n{traceback.format_exc()}")

    print(f"\n-- {tag}: one identical graph, different Leiden backends")
    for name, value in labels.items():
        describe(name, value)
    print()
    names = list(labels)
    for index, left in enumerate(names):
        for right in names[index + 1 :]:
            agreement(f"{left} vs {right}", labels[left], labels[right])


def main() -> None:
    adata = load_imc()
    print(f"imc: {adata.shape}, {adata.obs[IMC_CLUSTER_KEY].cat.categories.size} cell types")

    banner("1. End-to-end divergence, as the suite measures it")
    for flavor, key, kwargs in (
        (
            "neighborhood",
            f"nhood_niche_res={RESOLUTION}",
            {
                "groups": IMC_CLUSTER_KEY,
                "n_neighbors": N_NEIGHBORS,
                "scale": True,
                "abs_nhood": False,
                "distance": 1,
            },
        ),
        ("utag", f"utag_niche_res={RESOLUTION}", {"n_neighbors": N_NEIGHBORS}),
    ):
        cpu = sq.gr.calculate_niche(adata, flavor=flavor, resolutions=RESOLUTION, inplace=False, **kwargs)
        gpu = rsc.gr.calculate_niche(adata, flavor=flavor, resolutions=RESOLUTION, random_state=0, copy=True, **kwargs)
        print(f"\n-- {flavor}")
        describe("cpu", cpu.obs[key])
        describe("gpu", gpu.obs[key])
        agreement("cpu vs gpu", cpu.obs[key].astype(str), gpu.obs[key].astype(str))

    banner("2. Are the UTAG features and their PCA numerically equivalent?")
    # Per-component correlation is only meaningful where eigenvalues are separated. Where
    # they are nearly tied the individual directions are not identified, so the subspace
    # alignment is reported alongside.
    features = utag_feature_matrix(adata)
    cpu = AnnData(features.copy())
    sc.tl.pca(cpu)
    gpu = AnnData(features.copy())
    rsc.get.anndata_to_GPU(gpu)
    rsc.pp.pca(gpu)
    rsc.get.anndata_to_CPU(gpu)
    left = np.asarray(cpu.obsm["X_pca"])
    right = np.asarray(rsc.get.X_to_CPU(gpu.obsm["X_pca"]))
    variance = np.asarray(cpu.uns["pca"]["variance_ratio"])

    print(f"utag features: {features.shape}   pca: {left.shape} vs {right.shape}")
    print(f"{'comp':>4} {'|r|':>9} {'cpu var%':>9} {'cumulative%':>12}")
    cumulative = 0.0
    for index in range(min(left.shape[1], right.shape[1])):
        correlation = abs(np.corrcoef(left[:, index], right[:, index])[0, 1])
        share = float(variance[index]) * 100 if index < variance.size else float("nan")
        cumulative += share
        flag = "  <-- unidentified direction" if correlation < 0.99 else ""
        print(f"{index:>4} {correlation:>9.6f} {share:>9.4f} {cumulative:>12.4f}{flag}")

    for width in (2, 5, 10, 15, left.shape[1]):
        width = min(width, left.shape[1], right.shape[1])
        left_basis = np.linalg.qr(left[:, :width])[0]
        right_basis = np.linalg.qr(right[:, :width])[0]
        singular = np.linalg.svd(left_basis.T @ right_basis, compute_uv=False)
        print(f"first {width:2d} components: subspace alignment min={singular.min():.6f} (1.0 = identical span)")

    banner("3. Does kNN graph construction differ on identical features?")
    profile = scaled_neighborhood_profile(adata)
    print(f"scaled neighborhood profile: {profile.shape}")
    # Duplicate rows are the root cause of the disagreement: they create exact distance
    # ties, and with ties any of several neighbor sets is equally correct.
    distinct = np.unique(profile, axis=0).shape[0]
    print(f"  distinct profiles: {distinct} of {profile.shape[0]} -> {profile.shape[0] - distinct} duplicates")
    accuracy_against_exact(profile)
    cpu_graph = AnnData(profile.copy())
    sc.pp.neighbors(cpu_graph, n_neighbors=N_NEIGHBORS, use_rep="X", random_state=0)
    gpu_graph = AnnData(profile.copy())
    rsc.get.anndata_to_GPU(gpu_graph)
    rsc.pp.neighbors(gpu_graph, n_neighbors=N_NEIGHBORS, use_rep="X", algorithm="brute")
    rsc.get.anndata_to_CPU(gpu_graph)
    print(
        f"distances      jaccard cpu vs gpu = {graph_jaccard(cpu_graph.obsp['distances'], gpu_graph.obsp['distances']):.6f}"
    )
    print(
        "connectivities jaccard cpu vs gpu = "
        f"{graph_jaccard(cpu_graph.obsp['connectivities'], gpu_graph.obsp['connectivities']):.6f}"
    )

    # Rule out approximate search as the explanation before blaming tie-breaking.
    gpu_default = AnnData(profile.copy())
    rsc.get.anndata_to_GPU(gpu_default)
    rsc.pp.neighbors(gpu_default, n_neighbors=N_NEIGHBORS, use_rep="X")
    rsc.get.anndata_to_CPU(gpu_default)
    print(
        "gpu approximate vs gpu brute       = "
        f"{graph_jaccard(gpu_default.obsp['distances'], gpu_graph.obsp['distances']):.6f}"
    )

    banner("4. Leiden backends on identical graphs")
    compare_leiden_backends(profile, "neighborhood profile (scaled)")
    pca_input = AnnData(features.copy())
    sc.tl.pca(pca_input)
    compare_leiden_backends(pca_input.obsm["X_pca"], "utag PCA")


if __name__ == "__main__":
    main()
