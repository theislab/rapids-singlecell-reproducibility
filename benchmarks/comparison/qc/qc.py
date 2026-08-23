import numpy as np
import rapids_singlecell as rsc
import scanpy as sc

adata_sc = sc.datasets.pbmc3k()
adata_sc.var["mt"] = adata_sc.var_names.str.startswith("MT-")
adata_sc.var["ribo"] = adata_sc.var_names.str.startswith(("RPS", "RPL"))
adata_sc.var["hb"] = adata_sc.var_names.str.contains("^HB[^(P)]")
adata_rsc = adata_sc.copy()
rsc.get.anndata_to_GPU(adata_rsc)

sc.pp.calculate_qc_metrics(adata_sc, qc_vars=["mt", "ribo", "hb"], log1p=True, percent_top=False, inplace=True)
rsc.pp.calculate_qc_metrics(adata_rsc, qc_vars=["mt", "ribo", "hb"], log1p=True)

for key in adata_sc.obs.columns:
    np.testing.assert_allclose(adata_sc.obs[key], adata_rsc.obs[key], rtol=1e-6, atol=1e-6)
