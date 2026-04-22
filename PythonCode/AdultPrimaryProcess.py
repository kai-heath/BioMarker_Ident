# Core scverse libraries
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — renders to file only

import anndata as ad

# Data retrieval
import scanpy as sc
import scanpy.external as sce

# Reproducibility
import scvi

filePath = "RawData/HCA_Adult_Primary_Tissue/Global_raw.h5ad"
adata = sc.read_h5ad(filePath)

# paper using this data only used these donors, idk if their reasoning translates to my 
# usage, but it shrinks the size so...


# 3. Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

#sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
sc.tl.leiden(adata, resolution = 0.1, flavor="igraph", n_iterations=2, random_state=42)

sc.pl.umap(adata, size=4, show=False, save="_umapAPT.png")

adata = adata[adata.obs["leiden"].isin(["5", "6", "2"])].copy()
adata = adata[adata.obs["donor"].isin(["D2", "D3", "D4"])]
adata.write_zarr("ProcessedData/APTCardioScanpy.zarr")
