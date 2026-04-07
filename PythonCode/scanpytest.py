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

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)


sc.settings.set_figure_params(dpi=50, facecolor="white")

samples = {
    "Run1": "RawData/GSE202398_(D0-D30Diff)/GSM6118768_Galdos_Seq_Run1_filtered_feature_bc_matrix.h5",
    "Run2": "RawData/GSE202398_(D0-D30Diff)/GSM6118769_Galdos_Seq_Run2_filtered_feature_bc_matrix.h5",
    "Run3": "RawData/GSE202398_(D0-D30Diff)/GSM6118770_Galdos_Seq_Run3_filtered_feature_bc_matrix.h5",
}

adatas = {}
hto_cols = set()

for sample_id, filename in samples.items():
    # Load with gex_only=False to get both GEX and Antibody Capture
    sample_adata = sc.read_10x_h5(filename, gex_only=False)
    sample_adata.var_names_make_unique()

    # Split into GEX and Antibody Capture
    is_hto = sample_adata.var["feature_types"] == "Antibody Capture"
    gex_adata = sample_adata[:, ~is_hto].copy()
    hto_adata = sample_adata[:, is_hto].copy()

    # Collect HTO column names
    hto_cols.update(set(hto_adata.var_names))

    # Move HTO counts to obs in the gex_adata
    for i, hto_name in enumerate(hto_adata.var_names):
        # Convert to dense and flatten to match obs column expectations
        gex_adata.obs[hto_name] = hto_adata.X[:, i].toarray().flatten()

    adatas[sample_id] = gex_adata

hto_cols = list(hto_cols)

adata = ad.concat(adatas, label="sample", join="outer")
adata.obs_names_make_unique()

# Fill NaNs in HTO columns (from outer join) with 0
adata.obs[hto_cols] = adata.obs[hto_cols].fillna(0)

print(f"Total cells: {adata.n_obs}")
print(adata.obs["sample"].value_counts())

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

# Now QC metrics will only count genes
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    show=False,
    save="_qc_violin.png",
)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False, save="_qc_scatter.png")

sce.pp.hashsolo(adata, list(hto_cols))

adata = adata[~adata.obs['Classification'].isin(["Doublet", "Negative"]), :].copy()
adata.obs["CellLine"] = adata.obs["Classification"].str.extract(r"^([^_]+)")
adata.obs["Day"] = adata.obs["Classification"].str.extract(r"DAY(\d+)").astype(int).astype("category")

sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, min_counts=1000)

sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
#adata = adata[adata.obs.pct_counts_ribo < 20, :].copy()

# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(adata, show=False, save="_hvg.png")

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False, save="_pca_variance.png")

sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
    show=False,
    save="_pca.png",
)

sc.pp.neighbors(adata)

sc.tl.umap(adata, random_state=42)
sc.pl.umap(
    adata,
    color="Day",
    size=4,
    show=False,
    save="_umap_Day.png",
)

# Using the igraph implementation and a fixed number of iterations can be significantly faster,
# especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, random_state=42)

rcParams['figure.figsize'] = (6, 6)
sc.pl.umap(adata, color=["leiden"], show=False, save="_umap_leiden.png")

sc.pl.umap(
    adata,
    color=["CellLine", "Day"],
    wspace=0.5,
    size=3,
    show=False,
    save="_umap_leiden_doublets.png",
).astype(str).astype("category")

adata.write_h5ad("processed_adata.h5ad")

