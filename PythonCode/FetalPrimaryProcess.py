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

from pathlib import Path

directory_path = Path("RawData/Fetal_Primary/1_Chromium_cellranger_data_SC")
files = [f.name for f in directory_path.iterdir() if f.is_file()]
adatas = {}
for file in files:
    path = "RawData/Fetal_Primary/1_Chromium_cellranger_data_SC/" + file
    tmp_adata = sc.read_10x_h5(path)
    tmp_adata.var_names_make_unique()
    adatas[file] = tmp_adata

adata = ad.concat(adatas, label="sample_batch")
