import palantir
import scanpy as sc
import pandas as pd
import os
import anndata as ad
# Plotting
import matplotlib
import matplotlib.pyplot as plt

# warnings
import warnings
from numba.core.errors import NumbaDeprecationWarning

warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)

adata = ad.read_zarr("ProcessedData/mergedScanpy.zarr")

sc.pp.normalize_total(adata)
palantir.preprocess.log_transform(adata)


dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5, pca_key = "X_scVI")
ms_data = palantir.utils.determine_multiscale_space(adata)

sc.pl.embedding(
    adata,
    basis="umap",
    frameon=False,
)



terminal_states = pd.Series(
    ["Ventrical", "Atrial"],
    index=["HCAHeart7702880_AGGTCATTCATTCACT-1", "HCAHeart7833852_AGCTTGATCAACACCA"],
)

palantir.plot.highlight_cells_on_umap(adata, terminal_states)



imputed_X = palantir.utils.run_magic_imputation(adata, n_jobs=2)

sc.pl.embedding(
    adata,
    basis="umap",
    layer="MAGIC_imputed_data",
    color=["ENSG00000106631"],
    frameon=False,
)

palantir.plot.plot_diffusion_components(adata)

start_cell = "AAACCCAAGAGCATAT-1"

pr_res = palantir.core.run_palantir(
    adata, start_cell, num_waypoints=500, terminal_states=terminal_states,
)

palantir.plot.plot_palantir_results(adata, s=3)
