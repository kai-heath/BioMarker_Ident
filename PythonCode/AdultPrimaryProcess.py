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

adata = adata[adata.obs["donor"].isin(["D2", "D3", "D4"])].copy()

# 3. Normalize
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

#sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, random_state=42)

sc.pl.umap(adata, size=4, color = "scANVI_predictions", show=False, save="_umapAPT.png")

adata = adata[adata.obs["scANVI_predictions"].isin(["Ventricular Cardiomyocyte", "Atrial Cardiomyocyte"])].copy()

adata.write_zarr("ProcessedData/APTCardioScanpy.zarr")
