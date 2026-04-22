# Core scverse libraries
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — renders to file only

import anndata as ad
import numpy as np

# Data retrieval
import scanpy as sc
import scanpy.external as sce

# Reproducibility
import scvi
import torch

# Running soupx from R in python
import logging

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

import lamindb as ln
import numpy as np
import scanpy as sc
import seaborn as sns
import pandas as pd
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
from scipy.sparse import csc_matrix
from scipy.stats import median_abs_deviation
rcb.logger.setLevel(logging.ERROR)
scvi.settings.seed = 0

def DoSoup(adata, rawFile):
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp, target_sum=1e4)
    sc.pp.log1p(adata_pp)

    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(
    adata_pp, key_added="soupx_groups", flavor="igraph", n_iterations=2, directed=False
    )

    # Preprocess variables for SoupX
    adata.obs["soupx_groups"] = adata_pp.obs["soupx_groups"]
    del adata_pp

    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T

    adata_raw = sc.read_10x_h5(rawFile)
    adata_raw.var_names_make_unique()
    shared_genes = adata.var_names.intersection(adata_raw.var_names)
    adata = adata[:, shared_genes].copy()
    adata_raw = adata_raw[:, shared_genes].copy()
    # ensure same order
    adata_raw = adata_raw[:, adata.var_names].copy()


    genes_raw = adata_raw.var_names
    cells_raw = adata_raw.obs_names

    data_tod = adata_raw.X.T
    del adata_raw

    data_csc = data.tocsc()
    data_tod_csc = data_tod.tocsc()

    # Extract sparse components and cast to correct types
    x = data_csc.data.astype(np.float64)
    i = data_csc.indices.astype(np.int32)
    p = data_csc.indptr.astype(np.int32)
    dims = np.array(data_csc.shape, dtype=np.int32)

    x_tod = data_tod_csc.data.astype(np.float64)
    i_tod = data_tod_csc.indices.astype(np.int32)
    p_tod = data_tod_csc.indptr.astype(np.int32)
    dims_tod = np.array(data_tod_csc.shape, dtype=np.int32)

    with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        ro.globalenv["x"] = x
        ro.globalenv["i"] = i
        ro.globalenv["p"] = p
        ro.globalenv["dims"] = dims

        ro.globalenv["x_tod"] = x_tod
        ro.globalenv["i_tod"] = i_tod
        ro.globalenv["p_tod"] = p_tod
        ro.globalenv["dims_tod"] = dims_tod

        ro.globalenv["genes"] = np.array(genes)
        ro.globalenv["genes_raw"] = np.array(genes_raw)
        ro.globalenv["cells"] = np.array(cells)
        ro.globalenv["cells_raw"] = np.array(cells_raw)
        ro.globalenv["soupx_groups"] = adata.obs["soupx_groups"].to_numpy()

    ro.r ('''
    library(Matrix)
    library(SoupX)

    # Manually coerce types to avoid "array" class errors
    x <- as.numeric(x)
    i <- as.integer(i)
    p <- as.integer(p)
    dims <- as.integer(dims)

    x_tod <- as.numeric(x_tod)
    i_tod <- as.integer(i_tod)
    p_tod <- as.integer(p_tod)
    dims_tod <- as.integer(dims_tod)

    # Reconstruct sparse matrices
    data <- new("dgCMatrix",
            Dim = dims,
            x = x,
            i = i,
            p = p)

    data_tod <- new("dgCMatrix",
                Dim = dims_tod,
                x = x_tod,
                i = i_tod,
                p = p_tod)

    # Assign row and column names
    rownames(data) <- genes
    colnames(data) <- cells
    rownames(data_tod) <- genes_raw
    colnames(data_tod) <- cells_raw

    # SoupX pipeline
    sc = SoupChannel(data_tod, data, calcSoupProfile = TRUE)
    sc = setClusters(sc, soupx_groups)
    sc = autoEstCont(sc, doPlot = FALSE)
    out = adjustCounts(sc, roundToInt = TRUE)
    ''')
    with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        out_py = ro.conversion.rpy2py(ro.globalenv["out"])

    x = np.array(out_py.slots["x"])
    i = np.array(out_py.slots["i"])
    p = np.array(out_py.slots["p"])
    shape = tuple(out_py.slots["Dim"])

    out_matrix = csc_matrix((x, i, p), shape=shape)

    adata.layers["counts"] = adata.X.copy()
    adata.layers["soupX_counts"] = out_matrix.T
    adata.X = adata.layers["soupX_counts"]

    print(f"Total number of genes: {adata.n_vars}")

    # Min 20 cells - filters out 0 count genes
    sc.pp.filter_genes(adata, min_cells=20)
    print(f"Number of genes after cell filter: {adata.n_vars}")
    return adata

def is_outlier(adata, metric: str, nmads: int, upper_only: bool = False):
    M = adata.obs[metric]
    outlier = pd.Series(False, index=adata.obs_names, dtype=bool)

    for celltype in adata.obs["Classification"].unique():
        isCelltype = adata.obs["Classification"] == celltype
        M_group = M[isCelltype]

        med = np.median(M_group)
        mad = median_abs_deviation(M_group)

        if upper_only:
            group_outlier = M_group > med + nmads * mad
        else:
            group_outlier = (M_group < med - nmads * mad) | (M_group > med + nmads * mad)

        outlier[isCelltype] = group_outlier

    return outlier


Galdos_Hipsc = [
        "RawData/GSE202398_(D0-D30Diff)/GSM6118768_Galdos_Seq_Run1_filtered_feature_bc_matrix.h5",
        "RawData/GSE202398_(D0-D30Diff)/GSM6118769_Galdos_Seq_Run2_filtered_feature_bc_matrix.h5",
        "RawData/GSE202398_(D0-D30Diff)/GSM6118770_Galdos_Seq_Run3_filtered_feature_bc_matrix.h5",
    ]
rawFiles = [

    "RawData/GSE202398_(D0-D30Diff)/GSM6118768_Galdos_Seq_Run1_raw_feature_bc_matrix.h5",
    "RawData/GSE202398_(D0-D30Diff)/GSM6118769_Galdos_Seq_Run2_raw_feature_bc_matrix.h5",
    "RawData/GSE202398_(D0-D30Diff)/GSM6118770_Galdos_Seq_Run3_raw_feature_bc_matrix.h5",
]
HTOnames = [['WTC_DAY1', 'WTC_DAY2', 'WTC_DAY3', 'WTC_DAY4', 'WTC_DAY5', 'WTC_DAY6', 'LMNA_DAY1', 'LMNA_DAY2', 'LMNA_DAY3', 'LMNA_DAY4', 'LMNA_DAY5', 'LMNA_DAY6'],
["WTC_DAY0", "WTC_DAY7", "WTC_DAY11", "LMNA_DAY0", "LMNA_DAY7", "LMNA_DAY11"],
["WTC_DAY13", "WTC_DAY15", "WTC_DAY30", "LMNA_DAY13", "LMNA_DAY15", "LMNA_DAY30"]
]

GaldosAdatas = []
for i, filePath in enumerate(Galdos_Hipsc):
    # get both Gene Expression (gex) and Antibody Capture
    sample_adata = sc.read_10x_h5(filePath, gex_only=False)
    sample_adata.var_names_make_unique()

    # Split into Gene Expression (GEX) and Antibody Capture
    is_hto = sample_adata.var["feature_types"] == "Antibody Capture"
    gex_adata = sample_adata[:, ~is_hto].copy()
    hto_adata = sample_adata[:, is_hto].copy()
    del sample_adata
    
    # Move HTO counts to obs in the gex_adata
    for j, hto_name in enumerate(hto_adata.var_names):
        # Convert to dense and flatten to match obs column expectations
        gex_adata.obs[hto_name] = hto_adata.X[:, j].toarray().flatten()
    gex_adata.layers["counts"] = gex_adata.X.copy()
    sce.pp.hashsolo(gex_adata, list(HTOnames[i]))
    gex_adata = gex_adata[~gex_adata.obs['Classification'].isin(["Doublet", "Negative", "None"]), :].copy()
    
    # Auto QC
    # mitochondrial genes
    gex_adata.var["mt"] = gex_adata.var_names.str.startswith("MT-")
    # ribosomal genes
    gex_adata.var["ribo"] = gex_adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    gex_adata.var["hb"] = gex_adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
    sc.pp.calculate_qc_metrics(
        gex_adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    # defining what is an outlier based on median absolute deviation
    gex_adata.obs["outlier"] = (
    is_outlier(gex_adata, "log1p_total_counts", 5)
    | is_outlier(gex_adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(gex_adata, "pct_counts_in_top_20_genes", 5)
    )
    gex_adata.obs["mt_outlier"] = is_outlier(gex_adata, "pct_counts_mt", 5, True) 
    
    print(f"Total number of cells: {gex_adata.n_obs}")
    gex_adata = gex_adata[(~gex_adata.obs.outlier) & (~gex_adata.obs.mt_outlier)].copy()

    print(f"Number of cells after filtering of low quality cells: {gex_adata.n_obs}")
    gex_adata.raw = gex_adata
    gex_adata = DoSoup(gex_adata, rawFiles[i])

    GaldosAdatas.append(gex_adata)
    del gex_adata, hto_adata
    

GaldosAdatas = ad.concat(GaldosAdatas, label="sample", join="outer")
GaldosAdatas.var_names_make_unique()
GaldosAdatas.obs_names_make_unique()
sc.pp.highly_variable_genes(
    GaldosAdatas,
    n_top_genes=3000,
    subset=True,
    layer="soupX_counts",
    flavor="seurat_v3",
    batch_key="sample",
)
GaldosAdatas.obs["CellLine"] = GaldosAdatas.obs["Classification"].str.extract(r"^([^_]+)")
GaldosAdatas.obs["Day"] = GaldosAdatas.obs["Classification"].str.extract(r"DAY(\d+)").astype(int).astype("category")
toRemove = GaldosAdatas.obs.columns[GaldosAdatas.obs.columns.str.contains(r"_DAY")].tolist()
GaldosAdatas.obs.drop(columns = toRemove, inplace = True)

#sc.pp.pca(GaldosAdatas)

torch.set_float32_matmul_precision('high')
scvi.model.SCVI.setup_anndata(
    GaldosAdatas,
    layer="counts",
    batch_key = "sample",
    categorical_covariate_keys=["CellLine"],
    continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
)
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
    save="_umap_DayGaldos.png",
)

# 1. Convert .obs and .var row indices to strings
GaldosAdatas.obs_names = GaldosAdatas.obs_names.astype(str)
GaldosAdatas.var_names = GaldosAdatas.var_names.astype(str)

# 2. Convert .obs and .var column headers to strings
GaldosAdatas.obs.columns = GaldosAdatas.obs.columns.astype(str)
GaldosAdatas.var.columns = GaldosAdatas.var.columns.astype(str)

GaldosAdatas.write_zarr("ProcessedData/GaldosScanpy.zarr")

