import anndata as ad
import numpy as np
import pandas as pd
import pooch
import rapids_singlecell as rsc
from scipy.stats import pearsonr

X_pca_file = pooch.retrieve(
    "https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/pbmc_3500_pcs.tsv.gz",
    known_hash="md5:27e319b3ddcc0c00d98e70aa8e677b10",
)
X_pca = pd.read_csv(X_pca_file, delimiter="\t")

X_pca_harmony_file = pooch.retrieve(
    "https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/pbmc_3500_pcs_harmonized.tsv.gz",
    known_hash="md5:a7c4ce4b98c390997c66d63d48e09221",
)
X_pca_harmony_ref = pd.read_csv(X_pca_harmony_file, delimiter="\t").values

meta_file = pooch.retrieve(
    "https://github.com/slowkow/harmonypy/raw/refs/heads/master/data/pbmc_3500_meta.tsv.gz",
    known_hash="md5:8c7ca20e926513da7cf0def1211baecb",
)
meta = pd.read_csv(meta_file, delimiter="\t")

adata = ad.AnnData(X=None, obs=meta, obsm={"X_pca": X_pca.values})

rsc.pp.harmony_integrate(adata, key="donor", max_iter_harmony=20)

# Compare against harmonypy reference embedding per PC
corr = np.array(
    [pearsonr(adata.obsm["X_pca_harmony"][:, i], X_pca_harmony_ref[:, i])[0] for i in range(X_pca_harmony_ref.shape[1])]
)
l2 = np.linalg.norm(adata.obsm["X_pca_harmony"] - X_pca_harmony_ref, axis=0) / np.linalg.norm(X_pca_harmony_ref, axis=0)
assert corr.min() > 0.95
assert l2.max() < 0.1
