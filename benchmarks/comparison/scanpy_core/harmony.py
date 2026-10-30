from __future__ import annotations

import anndata as ad
import harmonypy
import pandas as pd
import pooch
import rapids_singlecell as rsc

from benchmarks.comparison.arrays import capture
from benchmarks.comparison.report import write_report

METHOD = "scanpy_core_harmony"

pcs_path = pooch.retrieve(
    "https://github.com/slowkow/harmonypy/raw/3cb114c/data/pbmc_3500_pcs.tsv.gz",
    known_hash="md5:27e319b3ddcc0c00d98e70aa8e677b10",
)
meta_path = pooch.retrieve(
    "https://github.com/slowkow/harmonypy/raw/3cb114c/data/pbmc_3500_meta.tsv.gz",
    known_hash="md5:8c7ca20e926513da7cf0def1211baecb",
)
pcs = pd.read_csv(pcs_path, delimiter="\t").to_numpy()
meta = pd.read_csv(meta_path, delimiter="\t")

reference = harmonypy.run_harmony(pcs, meta, "donor", max_iter_harmony=20, random_state=0).Z_corr
candidate = ad.AnnData(obs=meta, obsm={"X_pca": pcs})
rsc.pp.harmony_integrate(candidate, "donor", max_iter_harmony=20, random_state=0)

capture(METHOD, "harmony", reference=reference, candidate=candidate.obsm["X_pca_harmony"])
write_report(METHOD, "Harmonypy PBMC 3,500-cell donor benchmark", [], reference_package="harmonypy", shape=pcs.shape)
