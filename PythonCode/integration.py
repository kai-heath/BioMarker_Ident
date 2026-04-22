# Core scverse libraries
from mlx.core.random import categorical
from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt

import anndata as ad

# Data retrieval
import scanpy as sc
import scanpy.external as sce

# Reproducibility
import scvi
import torch

GaldosPath = "ProcessedData/GaldosScanpy.zarr"
BioreactorPath = "ProcessedData/BioReactorScanpy.zarr"


GaldosAdata = ad.read_zarr(GaldosPath)
BioreactorAdata = ad.read_zarr(BioreactorPath)

GaldosAdata.obs["source"] = "Galdos"
BioreactorAdata.obs["source"] = "Bioreactor"

BioreactorAdata.obs["Day"] = 15
BioreactorAdata.obs["Day"] = BioreactorAdata.obs["Day"].astype(int).astype("category")
BioreactorAdata.obs[['CellLine', 'substrate']] = BioreactorAdata.obs['orig.ident'].str.split('_', n=1, expand=True)

GaldosAdata.obs["substrate"] = "2D"


mergedAdata = ad.concat([GaldosAdata, BioreactorAdata], join = "inner")

del(GaldosAdata)
del(BioreactorAdata)

sc.pp.highly_variable_genes(
    mergedAdata,
    n_top_genes=3000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="source",
)
torch.set_float32_matmul_precision('medium')
scvi.model.SCVI.setup_anndata(
    mergedAdata,
    layer="counts",
    batch_key = "source",
    categorical_covariate_keys=["CellLine"],
    continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
)
model = scvi.model.SCVI(mergedAdata)
model.train()
mergedAdata.obsm["X_scVI"] = model.get_latent_representation()
mergedAdata.obsm["X_normalized_scVI"] = model.get_normalized_expression()  

sc.pp.neighbors(mergedAdata
    , use_rep="X_scVI"
    )
sc.tl.umap(mergedAdata)
sc.tl.leiden(mergedAdata)

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

p1 = sc.pl.umap(
    mergedAdata[(mergedAdata.obs["source"] == "Galdos").copy() & (mergedAdata.obs["Day"] == 15).copy()],
    size=4,
    ax=axs[0],
    color = "Day"
)
p2 = sc.pl.umap(
    mergedAdata[(mergedAdata.obs["source"] == "Bioreactor")],
    size=4,
    ax=axs[1],
    color = "substrate",
)


plt.tight_layout()
plt.show()
plt.savefig("figures/BioreactorGaldosComparison.png")

sc.pl.umap(
    mergedAdata,
    size = 4,
    color = "Day"
)

mergedAdata.write_zarr("ProcessedData/GaldosBioReactor.zarr")
mergedAdata = ad.read_zarr("ProcessedData/GaldosBioReactor.zarr")

APTPath = "ProcessedData/APTCardioScanpy.zarr"

ATPAdata = ad.read_zarr(APTPath)

ATPAdata.obs["Day"] = 100
ATPAdata.obs["Day"] = ATPAdata.obs["Day"].astype("category")
ATPAdata.obs["source"] = "PrimaryAdult"
ATPAdata.obs["CellType"] = ATPAdata.obs["scANVI_predictions"]


mergedAdata = ad.concat([mergedAdata, ATPAdata], join = "outer")
del ATPAdata

scvi.model.SCVI.setup_anndata(
    mergedAdata,
    layer="counts",
    batch_key = "source",
    #categorical_covariate_keys=["CellLine"],
    #continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
)
model = scvi.model.SCVI(mergedAdata)
model.train()
mergedAdata.obsm["X_scVI"] = model.get_latent_representation()
mergedAdata.obsm["X_normalized_scVI"] = model.get_normalized_expression()  

sc.pp.neighbors(mergedAdata,
    use_rep="X_scVI"
    )
sc.tl.umap(mergedAdata)
sc.tl.leiden(mergedAdata)

sc.pl.umap(
    mergedAdata,
    size = 4,
    color = "Day"
)
