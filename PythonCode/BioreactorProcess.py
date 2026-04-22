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
adata = adata[adata.obs["orig.ident"] != "AK2087_organoid", :].copy()

adata.obs[["cell_line", "condition"]] = adata.obs["orig.ident"].str.split("_", expand=True)

adata.layers["counts"] = adata.X.astype(int).copy()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False, save="_qc_scatter.png")

# 1. Save raw counts before any normalization
adata.layers["counts"] = adata.X.copy()

# 3. Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# 4. HVGs → scale → PCA → neighbors → UMAP
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
adataHVG = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adataHVG, max_value=10)
sc.pp.pca(adataHVG, n_comps=30)
sc.pp.neighbors(adataHVG, n_pcs=20)
sc.tl.umap(adataHVG, n_components=20)

sc.pl.umap(adataHVG, color="orig.ident", size=4, show=False, save="_umapBioReactor.png")

adata.write_zarr("ProcessedData/BioReactorScanpy.zarr")