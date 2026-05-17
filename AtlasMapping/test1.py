import scanpy as sc
import anndata
anndata.read = anndata.read_h5ad
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import pybiomart
import pandas as pd
import mygene
import scipy.sparse as sp

GaldosPath = "ProcessedData/GaldosScanpyUnmerged.zarr"
sourcePath = "RawData/Adult_Foetal_Integrated_Hearts_VKS2022_cellxgene.h5ad"
modelPath = "AtlasMapping/Atlas/"

adata = anndata.read_zarr(GaldosPath)
sourceAdata = anndata.read(sourcePath)
sourceAdata = sourceAdata.raw.to_adata()
sourceAdata.var_names = sourceAdata.var['index'].astype(str)
del sourceAdata.var['index']

mg = mygene.MyGeneInfo()

# 1. Get the list of symbols from your adata
gene_symbols = adata.var_names.tolist()

# 2. Query MyGene.info
# Use species='human' (or '9606')
results = mg.querymany(gene_symbols, scopes='symbol', fields='ensembl.gene', species='human')

# 3. Create a mapping dictionary
# Handle the case where a symbol maps to multiple IDs or none
mapping = {}
for res in results:
    if 'ensembl' in res:
        # Ensembl can return a list or a single string
        ensembl_data = res['ensembl']
        if isinstance(ensembl_data, list):
            mapping[res['query']] = ensembl_data[0]['gene']
        else:
            mapping[res['query']] = ensembl_data['gene']

# 4. Update adata
adata.var['gene_symbols'] = adata.var_names # Keep symbols for plotting
adata.var_names = [mapping.get(s, s) for s in adata.var_names]
adata.var_names_make_unique()


# 1. Get the list of genes the model expects
# We pull this directly from the source model's metadata
expected_genes = sourceAdata.var_names

# 2. Create a placeholder matrix of zeros (Cells x Expected Genes)
new_X = sp.csr_matrix((adata.shape[0], len(expected_genes)))

# 3. Create the new AnnData object
adata_padded = anndata.AnnData(X=new_X, obs=adata.obs)
adata_padded.var_names = expected_genes

# 4. Fill in the data for the genes you actually have
common_genes = adata.var_names.intersection(expected_genes)
adata_padded[:, common_genes].X = adata[:, common_genes].X

# 5. Copy over the rest of the metadata (obsm, etc.)
for key in adata.obsm_keys():
    adata_padded.obsm[key] = adata.obsm[key]

# 6. Now use the padded object for loading
vae = sca.models.SCVI.load_query_data(
    adata_padded,
    modelPath,
    freeze_dropout=True,
)


early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

vae.train(max_epochs=500, early_stopping=early_stopping_kwargs)

hipsc_latent = sc.AnnData(vae.get_latent_representation())

hipsc_latent.obs = vae.adata.obs.copy() 

reference_latent = sc.AnnData(vae.get_latent_representation())



sc.pp.neighbors(hipsc_latent, n_neighbors=8)
sc.pp.pca(hipsc_latent)
sc.tl.umap(hipsc_latent, init_pos='pca', random_state=42)
sc.tl.leiden(hipsc_latent)

sc.pl.umap(hipsc_latent,
           color=['Day'],
           frameon=False,
           wspace=0.6,
           )
