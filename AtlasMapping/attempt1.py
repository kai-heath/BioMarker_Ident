import scanpy as sc
import anndata
anndata.read = anndata.read_h5ad
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np

dataPath = "RawData/Adult_Foetal_Integrated_Hearts_VKS2022_cellxgene.h5ad"

adata = anndata.read(dataPath)
adata_raw = adata.raw.to_adata()
adata_raw.var_names = adata_raw.var['index'].astype(str)
del adata_raw.var['index']

sca.models.SCVI.setup_anndata(adata_raw, batch_key="sample")

vae = sca.models.SCVI(
    adata_raw,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
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

ref_path = 'AtlasMapping/Atlas/'
vae.save(ref_path, overwrite=True)

