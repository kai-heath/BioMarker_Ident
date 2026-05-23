# Core scverse libraries
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — renders to file only

import anndata as ad

# Data retrieval
import scanpy as sc
import scanpy.external as sce

# Reproducibility
import numpy as np

filePath = "RawData/HCA_Adult_Primary_Tissue/Global_raw.h5ad"
adata = sc.read_h5ad(filePath)

# grab the labels and only get
adata = adata[adata.obs["scANVI_predictions"].isin(["Ventricular Cardiomyocyte", "Atrial Cardiomyocyte"])].copy()
adata = adata[adata.obs["donor"].isin(["D1", "D2", "D3"])].copy()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.subsample(adata, n_obs=10000)
batch_counts = adata.obs["batch"].value_counts()
keep_batches = batch_counts[batch_counts >= 1000].index
adata = adata[adata.obs["batch"].isin(keep_batches)].copy()


adata.layers["counts"] = adata.X.copy()

#sc.pp.highly_variable_genes(
#    adata,
#    layer = "counts",
#    batch_key="batch",
#    n_top_genes = 6000,
#    subset = True,
#    flavor = "seurat_v3"
#)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write_zarr("ProcessedData/APTCardioScanpy.zarr")
adata = ad.read_zarr("ProcessedData/APTCardioScanpy.zarr")
