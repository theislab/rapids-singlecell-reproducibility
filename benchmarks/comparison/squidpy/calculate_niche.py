from __future__ import annotations

import rapids_singlecell as rsc
import squidpy as sq
from _report import lower_bound, upper_bound, write_report
from _shared import IMC_CLUSTER_KEY, load_imc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

adata = load_imc()
metrics = []


def compare_labels(flavor, reference, candidate, key, minimum_agreement, *, agreement="", count=""):
    reference_labels = reference.obs[key].astype(str)
    candidate_labels = candidate.obs[key].astype(str)
    metrics.extend(
        [
            lower_bound(
                f"{flavor}.adjusted_rand_index",
                adjusted_rand_score(reference_labels, candidate_labels),
                minimum_agreement,
            )
            | {"diagnosis": agreement},
            lower_bound(
                f"{flavor}.normalized_mutual_information",
                normalized_mutual_info_score(reference_labels, candidate_labels),
                minimum_agreement,
            )
            | {"diagnosis": agreement},
            upper_bound(
                f"{flavor}.cluster_count_difference",
                abs(reference_labels.nunique() - candidate_labels.nunique()),
                1,
            )
            | {"diagnosis": count},
        ]
    )


# Both Leiden-based flavors fail, and both criteria are left in place. `niche_divergence_diagnostic.py`
# separates the features, the kNN graph and the Leiden backend behind these numbers.
LEIDEN_COUNT_DIAGNOSIS = (
    "`resolution` does not carry the same meaning across Leiden implementations: on identical input "
    "cuGraph found 34 clusters where leidenalg found 41, and Scanpy's own two backends already "
    "differ by 1."
)
NEIGHBORHOOD_DIAGNOSIS = (
    "The neighborhood profile holds only 755 distinct rows across 4668 cells, so 89.67% of cells "
    "have an exact distance tie at the k-th neighbour. Against exact float64 ground truth the GPU "
    "kNN is exact and `sc.pp.neighbors` is not, on 743 rows: the divergence is the CPU reference, "
    "not the GPU."
)
UTAG_DIAGNOSIS = (
    "Squidpy calls `sc.tl.leiden` without a flavor, so CPU and GPU use different Leiden backends. "
    "On identical input Scanpy's own leidenalg and igraph backends agree only at ARI 0.5041, below "
    "this threshold, so the criterion measures backend choice rather than correctness."
)


# Neighborhood-profile niches. Squidpy's CPU implementation currently uses Scanpy's
# default seed (0) for neighbors and Leiden, so the GPU run uses that seed explicitly.
reference = sq.gr.calculate_niche(
    adata,
    flavor="neighborhood",
    groups=IMC_CLUSTER_KEY,
    n_neighbors=15,
    resolutions=0.5,
    scale=True,
    abs_nhood=False,
    distance=1,
    inplace=False,
)
candidate = rsc.gr.calculate_niche(
    adata,
    flavor="neighborhood",
    groups=IMC_CLUSTER_KEY,
    n_neighbors=15,
    resolutions=0.5,
    scale=True,
    abs_nhood=False,
    distance=1,
    random_state=0,
    copy=True,
)
compare_labels(
    "neighborhood",
    reference,
    candidate,
    key="nhood_niche_res=0.5",
    minimum_agreement=0.90,
    agreement=NEIGHBORHOOD_DIAGNOSIS,
    count=LEIDEN_COUNT_DIAGNOSIS,
)

# UTAG niches include PCA, kNN, and Leiden, so label agreement is evaluated rather
# than requiring arbitrary cluster identifiers to match exactly.
reference = sq.gr.calculate_niche(
    adata,
    flavor="utag",
    n_neighbors=15,
    resolutions=0.5,
    inplace=False,
)
candidate = rsc.gr.calculate_niche(
    adata,
    flavor="utag",
    n_neighbors=15,
    resolutions=0.5,
    random_state=0,
    copy=True,
)
compare_labels(
    "utag",
    reference,
    candidate,
    key="utag_niche_res=0.5",
    minimum_agreement=0.85,
    agreement=UTAG_DIAGNOSIS,
    count=LEIDEN_COUNT_DIAGNOSIS,
)

# CellCharter finishes with independently implemented Gaussian-mixture clustering.
# The same seed and initialization policy are used on both implementations.
cellcharter_adata = adata.copy()
cellcharter_adata.X = sparse.csr_matrix(cellcharter_adata.X)
reference = sq.gr.calculate_niche(
    cellcharter_adata,
    flavor="cellcharter",
    distance=3,
    aggregation="mean",
    n_components=6,
    random_state=42,
    inplace=False,
)
candidate = rsc.gr.calculate_niche(
    cellcharter_adata,
    flavor="cellcharter",
    distance=3,
    aggregation="mean",
    n_components=6,
    gmm_init="random_from_data",
    random_state=42,
    copy=True,
)
compare_labels(
    "cellcharter",
    reference,
    candidate,
    key="cellcharter_niche",
    minimum_agreement=0.80,
)

write_report(
    method="calculate_niche",
    dataset="squidpy.datasets.imc",
    tier="stochastic",
    metrics=metrics,
)
