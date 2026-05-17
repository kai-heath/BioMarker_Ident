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
#sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
adataHVG = adata[:, adata.var["highly_variable"]].copy()
#sc.pp.scale(adataHVG, max_value=10)
sc.pp.pca(adataHVG, n_comps=30)
sc.pp.neighbors(adataHVG, n_pcs=20)
sc.tl.umap(adataHVG, n_components=20)

sc.pl.umap(adataHVG, color="orig.ident", size=4, show=False, save="_umapBioReactor.png")

adata.write_zarr("ProcessedData/BioReactorScanpy.zarr")
adata = ad.read_zarr("ProcessedData/BioReactorScanpy.zarr")

mg = mygene.MyGeneInfo()
geneNames = adata.var_names
results = mg.querymany(
    geneNames.tolist(),
    scopes='symbol',
    fields='entrezgene,ensembl.gene',
    species='human',
    as_dataframe=True,
    returnall=True  # returns {'out': df, 'dup': [...], 'missing': [...]}
)

df = results['out']

# For dups: just keep the first hit (usually fine for broad CM filtering)
df = df[~df.index.duplicated(keep='first')]

def get_ensembl(val):
    if isinstance(val, dict):
        return val.get('gene')
    elif isinstance(val, list):
        return val[0].get('gene')  # take first if multiple
    return None

df['ensembl_id'] = df['ensembl.gene'].apply(get_ensembl)

# Map onto adata.var
adata.var['ensembl_id'] = adata.var_names.map(df['ensembl.gene'])
adata.var['entrez_id'] = adata.var_names.map(df['entrezgene'])

# Sanity check
n_mapped = adata.var['entrez_id'].notna().sum()
print(f"Mapped: {n_mapped} / {adata.n_vars} genes")
adata.var['entrez_id']

adata = adata[:, adata.var['ensembl_id'].notna()].copy()
adata.var_names = adata.var['ensembl_id']
