# Core scverse libraries
from __future__ import annotations

import matplotlib

import anndata as ad

# Data retrieval
import scanpy as sc
import scanpy.external as sce

# Reproducibility
import scvi

import mygene

filePath = "RawData/GSE263372_(BioReator)/converted_data.h5ad"
adata = sc.read_h5ad(filePath)
adata.obs.rename(columns={"orig.ident" : "condition"}, inplace = True)
adata = adata[adata.obs["condition"] != "AK2087_organoid"]
adata.obs["condition"] = (adata.obs["condition"]).astype(str) + "_DAY15"
adata.obs["condition"] = adata.obs["condition"].astype("category")

adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 1. Save raw counts before any normalization
adata.layers["counts"] = adata.X.copy()

# 3. Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#sc.pp.highly_variable_genes(
#    adata,
#    layer = "counts",
#    batch_key="old.ident",
#    subset = True,
#    flavor = "seurat_v3"
#)

mg = mygene.MyGeneInfo()
geneNames = adata.var_names
results = mg.querymany(
geneNames.tolist(),
scopes='symbol',
fields='entrezgene,ensembl.gene',
species='human',
as_dataframe=True,
returnall=True 
)

df = results['out']

df = df[~df.index.duplicated(keep='first')]
adata.var['ensembl_id'] = adata.var_names.map(df['ensembl.gene'])
adata.var["gene_id"] = adata.var_names
del df
adata = adata[:, adata.var['ensembl_id'].notna()].copy()

adata.var_names = adata.var['ensembl_id']

adata.write_zarr("ProcessedData/BioReactorScanpy.zarr")
