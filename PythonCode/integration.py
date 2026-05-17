# Core scverse libraries
import pandas as pd
import regex
from scvi.model import SCANVI
from jax.numpy import log2
import anndata as ad
from anndata import read_zarr
from __future__ import annotations


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scvi.external import SysVI

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

BioreactorAdata.obs["source"] = "Bioreactor"
GaldosAdata.obs["source"] = "Galdos"

GaldosAdata.obs["substrate"] = "2D"
GaldosAdata.obs.rename(columns = {"Classification" : "condition"}, inplace = True)
BioreactorAdata.obs["batch"] = "Bioreactor"
BioreactorAdata.obs["batch"] = BioreactorAdata.obs["batch"].astype(str).astype("category")

mergedAdata = ad.concat([GaldosAdata, BioreactorAdata], join = "inner")

del(GaldosAdata)
del(BioreactorAdata)

#mergedAdata.obs["Day"] = mergedAdata.obs["Day"].astype(str).astype('category')

mergedAdata.write_zarr("ProcessedData/GaldosBioReactor.zarr")
mergedAdata = ad.read_zarr("ProcessedData/GaldosBioReactor.zarr")

APTPath = "ProcessedData/APTCardioScanpy.zarr"

ATPAdata = ad.read_zarr(APTPath)

ATPAdata.obs["source"] = "PrimaryAdult"
ATPAdata.obs.rename(columns = {"scANVI_predictions" : "condition"}, inplace = True)
ATPAdata.obs["system"] = "ATCM"
mergedAdata.obs["system"] = "hipsc"

mergedAdata = ad.concat([mergedAdata, ATPAdata], join = "inner")
del ATPAdata

mergedAdata.obs["Day"] = ""
for i in enumerate(mergedAdata.obs["condition"]):
    condition = i[1]
    i = i[0]
    Day = regex.search(".*(DAY.*)", condition)
    if Day is not None:
        mergedAdata.obs["Day"][i] = Day[1]
    else: 
        mergedAdata.obs["Day"][i] = "cardiomyocyte"

mergedAdata.write_zarr("ProcessedData/GBA.zarr")
mergedAdata = ad.read_zarr("ProcessedData/GBA.zarr")

mergedAdata.obs["Day"] = mergedAdata.obs["Day"].astype(str).astype("category")
mergedAdata.obs["batch"] = mergedAdata.obs["batch"].astype(str).astype("category")

sc.pp.highly_variable_genes(
    mergedAdata,
    n_top_genes = 1500,
    batch_key = "batch"
)

scvi.model.SCVI.setup_anndata(mergedAdata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(mergedAdata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

SCVI_LATENT_KEY = "X_scVI"
mergedAdata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
mergedAdata.layers["SCVI_normalized"] = model.get_normalized_expression()
sc.pp.neighbors(mergedAdata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(mergedAdata, resolution = 0.1)

sc.tl.umap(mergedAdata)
sc.pl.umap(
    mergedAdata,
    color=["leiden", "Day"],
    frameon=False,
    ncols=1,
)

mergedAdata.write_zarr("ProcessedData/mergedScanpy.zarr")
mergedAdata.layers["scvi_log"] = np.log1p(mergedAdata.layers["SCVI_normalized"])
sc.tl.rank_genes_groups(
    mergedAdata, 
    groupby='leiden', 
    groups=['0'], 
    reference='3', 
    method='wilcoxon',
    use_raw=False,
    layer = "scvi_log"
)

annot = sc.queries.biomart_annotations(
    "hsapiens", 
    ["ensembl_gene_id", "external_gene_name", "description"], 
    host="www.ensembl.org"
).set_index("ensembl_gene_id")
annot = annot.drop_duplicates(subset=['external_gene_name'], keep='first')

mergedAdata.var['gene_symbols'] = mergedAdata.var.index.map(annot['external_gene_name'])

sc.pl.rank_genes_groups(mergedAdata, n_genes=20, sharey=False, gene_symbols="gene_symbols")

result_df = sc.get.rank_genes_groups_df(mergedAdata, group="0")

true_markers = result_df[
    (abs(result_df['logfoldchanges']) > 1.5) & 
    (result_df['pvals_adj'] < 1e-5)
]

true_markers["abslfc"] = abs(true_markers["logfoldchanges"]).copy()
df_sorted = true_markers.sort_values(by='abslfc', ascending=False)


df_sorted.to_csv("diffMarkers.csv")

maturity_markers = ['TNNI1', 'TNNI3', 'MYH6', 'MYH7', 'PLN', 'RYR2']

mergedAdata.var_names = mergedAdata.var['gene_symbols']

sc.pl.violin(
    mergedAdata, 
    maturity_markers, 
    groupby='leiden', 
    layer='SCVI_normalized'
)


sc.tl.rank_genes_groups(
    mergedAdata, 
    groupby='leiden', 
    groups=['1', '3', '5'], 
    reference='0', 
    method='wilcoxon',
    use_raw=False,
    layer = "scvi_log"
)

result_df = sc.get.rank_genes_groups_df(mergedAdata, group="5")
true_markers = result_df[
    (abs(result_df['logfoldchanges']) > 1.5) & 
    (result_df['pvals_adj'] < 1e-5)
]
true_markers["abslfc"] = abs(true_markers["logfoldchanges"]).copy()
df_sorted = true_markers.sort_values(by='abslfc', ascending=False)


df_sorted.to_csv("VentCMs.csv")

result_df = sc.get.rank_genes_groups_df(mergedAdata, group="1")
true_markers = result_df[
    (abs(result_df['logfoldchanges']) > 1.5) & 
    (result_df['pvals_adj'] < 1e-5)
]
true_markers["abslfc"] = abs(true_markers["logfoldchanges"]).copy()
earlyCM = true_markers.sort_values(by='abslfc', ascending=False)


earlyCM.to_csv("earlyCMs.csv")

result_df = sc.get.rank_genes_groups_df(mergedAdata, group="3")
true_markers = result_df[
    (abs(result_df['logfoldchanges']) > 1.5) & 
    (result_df['pvals_adj'] < 1e-5)
]
true_markers["abslfc"] = abs(true_markers["logfoldchanges"]).copy()
ATRCMs = true_markers.sort_values(by='abslfc', ascending=False)


ATRCMs.to_csv("AtrCMs.csv")

earlyToATR = ATRCMs.merge(earlyCM, on= 'names')
earlyToATR = earlyToATR[~ earlyToATR['names'].str.startswith('LINC')]


opposite_signs_df = earlyToATR[earlyToATR['logfoldchanges_y'] * earlyToATR['logfoldchanges_x'] < 0].copy()


opposite_signs_df['combinedChange'] = opposite_signs_df['abslfc_x'] * opposite_signs_df['abslfc_y']

opposite_signs_df = opposite_signs_df.sort_values(by = 'combinedChange', ascending=False)

opposite_signs_df.to_csv("analysis.csv")

refinedAnalysis = pd.DataFrame()
refinedAnalysis['GeneID'] = opposite_signs_df['names']
refinedAnalysis['logfc_Native_CM'] = opposite_signs_df['logfoldchanges_x']
refinedAnalysis['logfc_Early_hipsc_CM'] = opposite_signs_df['logfoldchanges_y']
refinedAnalysis['description'] = refinedAnalysis['GeneID'].map(annot.set_index('external_gene_name')['description'])
refinedAnalysis.to_csv("analysis.csv")
