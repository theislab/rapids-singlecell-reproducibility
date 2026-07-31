from __future__ import annotations

import rapids_singlecell as rsc
import squidpy as sq
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from _report import lower_bound, upper_bound, write_report
from _shared import IMC_CLUSTER_KEY, load_imc


adata = load_imc()
metrics = []


def compare_labels(flavor, reference, candidate, key, minimum_agreement):
    reference_labels = reference.obs[key].astype(str)
    candidate_labels = candidate.obs[key].astype(str)
    metrics.extend(
        [
            lower_bound(
                f"{flavor}.adjusted_rand_index",
                adjusted_rand_score(reference_labels, candidate_labels),
                minimum_agreement,
            ),
            lower_bound(
                f"{flavor}.normalized_mutual_information",
                normalized_mutual_info_score(reference_labels, candidate_labels),
                minimum_agreement,
            ),
            upper_bound(
                f"{flavor}.cluster_count_difference",
                abs(reference_labels.nunique() - candidate_labels.nunique()),
                1,
            ),
        ]
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
)

# CellCharter finishes with independently implemented Gaussian-mixture clustering.
# The same seed and initialization policy are used on both implementations.
reference = sq.gr.calculate_niche(
    adata,
    flavor="cellcharter",
    distance=3,
    aggregation="mean",
    n_components=6,
    random_state=42,
    inplace=False,
)
candidate = rsc.gr.calculate_niche(
    adata,
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
