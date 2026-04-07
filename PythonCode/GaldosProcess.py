# Core scverse libraries
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — renders to file only

import anndata as ad

# Data retrieval
import scanpy as sc
import scanpy.external as sce
from matplotlib import rcParams

# Reproducibility
import scvi

from PythonCode.utils import ScanpyCleaner
import harmonypy as hm
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)


Galdos_Hipsc = [
        "RawData/GSE202398_(D0-D30Diff)/GSM6118768_Galdos_Seq_Run1_filtered_feature_bc_matrix.h5",
        "RawData/GSE202398_(D0-D30Diff)/GSM6118769_Galdos_Seq_Run2_filtered_feature_bc_matrix.h5",
        "RawData/GSE202398_(D0-D30Diff)/GSM6118770_Galdos_Seq_Run3_filtered_feature_bc_matrix.h5",
    ]
RawData = [
    "RawData/GSE202398_(D0-D30Diff)/GSM6118768_Galdos_Seq_Run1_raw_feature_bc_matrix.h5",
        "RawData/GSE202398_(D0-D30Diff)/GSM6118769_Galdos_Seq_Run2_raw_feature_bc_matrix.h5",
        "RawData/GSE202398_(D0-D30Diff)/GSM6118770_Galdos_Seq_Run3_raw_feature_bc_matrix.h5"
]
# sample 1; Galdoes et al
GaldosAdatas = []
HTOnames = [['WTC_DAY1', 'WTC_DAY2', 'WTC_DAY3', 'WTC_DAY4', 'WTC_DAY5', 'WTC_DAY6', 'LMNA_DAY1', 'LMNA_DAY2', 'LMNA_DAY3', 'LMNA_DAY4', 'LMNA_DAY5', 'LMNA_DAY6'],
["WTC_DAY0", "WTC_DAY7", "WTC_DAY11", "LMNA_DAY0", "LMNA_DAY7", "LMNA_DAY11"],
["WTC_DAY13", "WTC_DAY15", "WTC_DAY30", "LMNA_DAY13", "LMNA_DAY15", "LMNA_DAY30"]
]
for i, filePath in enumerate(Galdos_Hipsc):
    HTOname = HTOnames[i]
    # get both Gene Expression (gex) and Antibody Capture
    sample_adata = sc.read_10x_h5(filePath, gex_only=False)
    raw_adata = sc.read_10x_h5(RawData[i])
    sample_adata.var_names_make_unique()
    raw_adata.var_names_make_unique()

    # Split into GEX and Antibody Capture
    is_hto = sample_adata.var["feature_types"] == "Antibody Capture"
    gex_adata = sample_adata[:, ~is_hto].copy()
    hto_adata = sample_adata[:, is_hto].copy()
    del sample_adata

    scvi.external.SCAR.setup_anndata(gex_adata)
    raw_adata = raw_adata[:, gex_adata.var_names].copy()
    scvi.external.SCAR.get_ambient_profile(adata=gex_adata, raw_adata=raw_adata, prob=0.995, n_batch=3)   
    vae = scvi.external.SCAR(gex_adata, ambient_profile="ambient_profile")
    vae.train()
    gex_adata.obsm["X_scAR"] = vae.get_latent_representation()
    gex_adata.layers['denoised'] = vae.get_denoised_counts()

    # Move HTO counts to obs in the gex_adata
    for j, hto_name in enumerate(hto_adata.var_names):
        # Convert to dense and flatten to match obs column expectations
        gex_adata.obs[hto_name] = hto_adata.X[:, j].toarray().flatten()
    sc.pp.log1p(gex_adata)
    sce.pp.hashsolo(gex_adata, list(HTOname))
    gex_adata = gex_adata[~gex_adata.obs['Classification'].isin(["Doublet", "Negative", "None"]), :].copy()
    sc.pp.filter_cells(gex_adata, min_genes=500)
    sc.pp.filter_genes(gex_adata, min_cells=3)
    sc.pp.filter_cells(gex_adata, min_counts=1000)
    gex_adata.var["mt"] = gex_adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(gex_adata, qc_vars=["mt"], inplace=True)
    gex_adata = gex_adata[gex_adata.obs["pct_counts_mt"] < 20, :].copy()

    sc.pp.normalize_total(gex_adata)
    GaldosAdatas.append(gex_adata)
    del gex_adata, hto_adata
    

GaldosAdatas = ad.concat(GaldosAdatas, label="sample", join="outer")
sc.pp.highly_variable_genes(GaldosAdatas)
GaldosAdatas.obs["CellLine"] = GaldosAdatas.obs["Classification"].str.extract(r"^([^_]+)")
GaldosAdatas.obs["Day"] = GaldosAdatas.obs["Classification"].str.extract(r"DAY(\d+)").astype(int).astype("category")
toRemove = GaldosAdatas.obs.columns[GaldosAdatas.obs.columns.str.contains(r"_DAY")].tolist()
GaldosAdatas.obs.drop(columns = toRemove, inplace = True)

sc.pp.pca(GaldosAdatas)

scvi.model.SCVI.setup_anndata(GaldosAdatas, batch_key="sample")
model = scvi.model.SCVI(GaldosAdatas)
model.train()
GaldosAdatas.obsm["X_scVI"] = model.get_latent_representation()
GaldosAdatas.obsm["X_normalized_scVI"] = model.get_normalized_expression()  
sc.pp.neighbors(GaldosAdatas, use_rep="X_scVI")
sc.tl.umap(GaldosAdatas)
sc.tl.leiden(GaldosAdatas)
sc.pl.umap(
    GaldosAdatas,
    color="Day",
    size=4,
    show=False,
    save="_umap_DayTEST.png",
)
