# Core scverse libraries
from anndata import read_zarr
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

mergedAdata = ad.concat([GaldosAdata, BioreactorAdata], join = "outer")

del(GaldosAdata)
del(BioreactorAdata)

mergedAdata.write_zarr("ProcessedData/GaldosBioReactor.zarr")
mergedAdata = ad.read_zarr("ProcessedData/GaldosBioReactor.zarr")

APTPath = "ProcessedData/APTCardioScanpy.zarr"

ATPAdata = ad.read_zarr(APTPath)

ATPAdata.obs["Day"] = 100
ATPAdata.obs["Day"] = ATPAdata.obs["Day"].astype("category")
ATPAdata.obs["source"] = "PrimaryAdult"
ATPAdata.obs["CellType"] = ATPAdata.obs["scANVI_predictions"]
ATPAdata.obs["cellOrigin"] = "native"

mergedAdata.obs["cellOrigin"] = "hipsc"
mergedAdata = ad.concat([mergedAdata, ATPAdata], join = "outer")
del ATPAdata

#sc.pp.normalize_total(mergedAdata, target_sum=1e4)
#sc.pp.log1p(mergedAdata)
sc.pp.highly_variable_genes(
    mergedAdata,
    n_top_genes=3000,
    subset=True,
    batch_key="source",
    layer="counts",
    flavor = "seurat_v3"
)

mergedAdata.obs["cellOrigin"] = mergedAdata.obs["cellOrigin"].astype("category")

torch.set_float32_matmul_precision('high')
scvi.model.SCVI.setup_anndata(
    mergedAdata,
    batch_key = "source",
    layer = "counts",
    #categorical_covariate_keys=["cellOrigin"],
    #continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
)
model = scvi.model.SCVI(mergedAdata)
model.train()
mergedAdata.obsm["X_scVI"] = model.get_latent_representation()
mergedAdata.obsm["X_normalized_scVI"] = model.get_normalized_expression()  

sc.pp.neighbors(mergedAdata,
    use_rep="X_scVI",
    #n_neighbors=10,
    metric="cosine"
    )
sc.tl.umap(
    mergedAdata,          # Increases whitespace between distinct groups
)

#sc.tl.leiden(mergedAdata, flavor='igraph', n_iterations=2)


sc.pl.umap(
    mergedAdata,
    size = 4,
    color = "substrate",
    save = "integratedsubstrate.png"
)

mergedAdata.write_zarr("ProcessedData/mergedScanpy.zarr")
mergedAdata = ad.read_zarr("ProcessedData/mergedScanpy.zarr")
