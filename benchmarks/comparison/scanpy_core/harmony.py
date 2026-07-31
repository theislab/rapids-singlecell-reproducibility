from __future__ import annotations

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
import pooch
import rapids_singlecell as rsc
from _report import lower, upper, write_report
from _shared import component_abs_correlations, max_abs

pcs_path = pooch.retrieve(
    "https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/pbmc_3500_pcs.tsv.gz",
    known_hash="md5:27e319b3ddcc0c00d98e70aa8e677b10",
)
meta_path = pooch.retrieve(
    "https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/pbmc_3500_meta.tsv.gz",
    known_hash="md5:8c7ca20e926513da7cf0def1211baecb",
)
pcs = pd.read_csv(pcs_path, delimiter="\t").to_numpy()
meta = pd.read_csv(meta_path, delimiter="\t")
reference = ad.AnnData(X=None, obs=meta.copy(), obsm={"X_pca": pcs.copy()})
candidate = reference.copy()

reference_result = harmonypy.run_harmony(
    np.asarray(reference.obsm["X_pca"], dtype=np.float64),
    reference.obs,
    "donor",
    max_iter_harmony=20,
    random_state=0,
    device="cpu",
)
reference_embedding = np.asarray(reference_result.Z_corr)
if reference_embedding.shape != reference.obsm["X_pca"].shape:
    reference_embedding = reference_embedding.T
if reference_embedding.shape != reference.obsm["X_pca"].shape:
    raise RuntimeError(
        f"Unexpected Harmonypy output shape {reference_embedding.shape}; expected {reference.obsm['X_pca'].shape}"
    )
reference.obsm["X_pca_harmony"] = reference_embedding

rsc.pp.harmony_integrate(
    candidate,
    key="donor",
    basis="X_pca",
    adjusted_basis="X_pca_harmony",
    flavor="harmony1",
    max_iter_harmony=20,
    max_iter_clustering=20,
    random_state=0,
)
rsc.get.anndata_to_CPU(candidate)

component_correlations = component_abs_correlations(reference.obsm["X_pca_harmony"], candidate.obsm["X_pca_harmony"])
metrics = [
    lower("harmony.minimum_component_abs_correlation", component_correlations.min(), 0.95),
    lower("harmony.mean_component_abs_correlation", component_correlations.mean(), 0.98),
    upper(
        "harmony.standard_deviation_max_abs_error",
        max_abs(reference.obsm["X_pca_harmony"].std(axis=0), candidate.obsm["X_pca_harmony"].std(axis=0)),
        0.1,
    ),
]

write_report(
    "scanpy_core_harmony",
    "Harmonypy PBMC 3,500-cell donor benchmark",
    "iterative",
    metrics,
    reference_package="harmonypy",
    packages=("harmonypy",),
)
