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

filePath = "RawData/GSE263372_(BioReator)/converted_data.h5ad"
adata = sc.read_h5ad(filePath)
adata

adata.layers["counts"] = adata.X.copy()

# Standard Scanpy preprocessing
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.raw = adata

