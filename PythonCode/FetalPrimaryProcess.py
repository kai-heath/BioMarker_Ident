# Core scverse libraries
from __future__ import annotations
import matplotlib
import anndata as ad
import pandas as pd# Data retrieval
import scanpy as sc
import numpy as np
import scanpy.external as sce
import h5py
from scipy.stats import median_abs_deviation

from pathlib import Path



def is_outlier(adata, metric: str, nmads: int, upper_only: bool = False):
    M = adata.obs[metric]
    outlier = pd.Series(False, index=adata.obs_names, dtype=bool)

    for celltype in adata.obs["Tissue"].unique():
        isCelltype = adata.obs["Tissue"] == celltype
        M_group = M[isCelltype]

        med = np.median(M_group)
        mad = median_abs_deviation(M_group)

        if upper_only:
            group_outlier = M_group > med + nmads * mad
        else:
            group_outlier = (M_group < med - nmads * mad) | (M_group > med + nmads * mad)

        outlier[isCelltype] = group_outlier

    return outlier





directory_path = Path("RawData/Fetal_Primary/1_Chromium_cellranger_data_SC")
files = [f.name for f in directory_path.iterdir() if f.is_file()][0:5]
adatas = {}
for file in files:
    path = "RawData/Fetal_Primary/1_Chromium_cellranger_data_SC/" + file
    print(path)
    i = h5py.File(path)['shoji']
    exp = pd.DataFrame(i["Expression"][()])
    cells = [c.decode('utf-8') for c in i["Cellid"][()]]
    genes = [g.decode('utf-8').split('.')[0] for g in i["Accession"][()]]
    geneNames = [h.decode('utf-8') for h in i["Gene"][()]]
    tissue = [t.decode('utf-8') for t in i["Tissue"][()]]
    doublet = list(i["DoubletFlag"][()])
    age = [t.astype(str) + "-fetal" for t in i["Age"][()]]
    adata = sc.AnnData(exp, obs = pd.DataFrame(index=cells), var = pd.DataFrame(index=genes))
    adata.var["Gene"] = geneNames
    adata.var_names_make_unique()
    adata.obs["Tissue"] = tissue
    adata.obs["DoubletFlag"] = doublet
    adata.obs["Day"] = age
    adata = adata[adata.obs["Tissue"] == "Heart"]

    adata = adata[adata.obs["DoubletFlag"] == np.False_]
    adata.var["mt"] = adata.var["Gene"].str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var["Gene"].str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var["Gene"].str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 5, True) 
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    adatas[file] = adata



adata = ad.concat(adatas, label="sample_batch")
sc.pp.filter_genes(adata, min_cells=3)
del adatas
adata.obs["cellOrigin"] = "fetal"
adata.obs["source"] = "fetal"
adata.layers["counts"] = adata.X.copy()

adata.write_zarr("ProcessedData/FTPScanpy.zarr")
