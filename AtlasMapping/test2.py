import scanpy as sc
import anndata
anndata.read = anndata.read_h5ad
import scarches as sca
import scvi
import scipy.sparse as sp
import mygene

GaldosPath = "ProcessedData/GaldosScanpyUnmerged.zarr"
sourcePath = "RawData/Adult_Foetal_Integrated_Hearts_VKS2022_cellxgene.h5ad"
modelPath = "AtlasMapping/Atlas/"

# --- Load data ---
adata = anndata.read_zarr(GaldosPath)
sourceAdata = anndata.read_h5ad(sourcePath)
sourceAdata = sourceAdata.raw.to_adata()
sourceAdata.var_names = sourceAdata.var['index'].astype(str)
del sourceAdata.var['index']

# --- Convert hiPSC gene symbols to Ensembl IDs ---
mg = mygene.MyGeneInfo()
results = mg.querymany(adata.var_names.tolist(), 
                       scopes='symbol', fields='ensembl.gene', species='human')
mapping = {}
for res in results:
    if 'ensembl' in res:
        ensembl_data = res['ensembl']
        mapping[res['query']] = (ensembl_data[0]['gene'] 
                                  if isinstance(ensembl_data, list) 
                                  else ensembl_data['gene'])

adata.var['gene_symbols'] = adata.var_names
adata.var_names = [mapping.get(s, s) for s in adata.var_names]
adata.var_names_make_unique()

# --- Check gene overlap before anything else ---
common_genes = adata.var_names.intersection(sourceAdata.var_names)
print(f"Reference genes: {sourceAdata.n_vars}")
print(f"Query genes: {adata.n_vars}")
print(f"Common genes: {len(common_genes)}")
# You want >3000 common genes — if much lower, the Ensembl conversion has issues

# --- Load reference model ---
vae_ref = sca.models.SCVI.load(modelPath, adata=sourceAdata)

# 1. Keep only genes that successfully converted to Ensembl IDs
ensembl_mask = adata.var_names.str.startswith('ENSG')
adata_ensembl = adata[:, ensembl_mask].copy()
print(f"After keeping only ENSG IDs: {adata_ensembl.n_vars} genes")

# 2. Check overlap with reference
common_genes = adata_ensembl.var_names.intersection(sourceAdata.var_names)
print(f"Common genes with reference: {len(common_genes)}")
# Want >3000 — if not, see note below

# 3. Subset to common genes only
adata_query = adata_ensembl[:, common_genes].copy()
print(f"Query adata ready: {adata_query.shape}")

# 4. Add required covariates matching reference setup
print(vae_ref.registry_['setup_args'])  # check what batch key was used
# Then add matching column:
adata_query.obs['sample'] = 'hiPSC_query'       # adjust key to match reference
adata_query.obs['suspension_type'] = 'cell'

# 5. Now load query data
vae_q = sca.models.SCVI.load_query_data(
    adata_query,
    reference_model=vae_ref,
    inplace_subset_query_vars=True,
)

# --- Train query model ---
vae_q.train(
    max_epochs=200,
    plan_kwargs={"weight_decay": 0.0},  # essential for scArches surgery
    early_stopping=True,
    early_stopping_patience=10,
)

# --- Get latent representations ---
reference_latent = sc.AnnData(vae_ref.get_latent_representation())
reference_latent.obs = vae_ref.adata.obs.copy()
reference_latent.obsm["X_scVI"] = reference_latent.X.copy()

hipsc_latent = sc.AnnData(vae_q.get_latent_representation())
hipsc_latent.obs = vae_q.adata.obs.copy()
hipsc_latent.obsm["X_scVI"] = hipsc_latent.X.copy()

# --- Build combined object for joint UMAP ---
reference_latent.obs['dataset'] = 'reference'
hipsc_latent.obs['dataset'] = 'query'

combined = sc.concat([reference_latent, hipsc_latent])
combined.obsm["X_scVI"] = np.vstack([
    reference_latent.obsm["X_scVI"],
    hipsc_latent.obsm["X_scVI"]
])

sc.pp.neighbors(combined, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(combined, init_pos='random', random_state=42)

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.umap(combined[combined.obs.dataset == 'reference'],
           color='cell_type', ax=axes[0], show=False, title='Reference cell types')

sc.pl.umap(combined, color='dataset',
           ax=axes[1], show=False, title='Reference vs Query')

sc.pl.umap(combined[combined.obs.dataset == 'query'],
           color='Day', ax=axes[2], show=False, title='hiPSC-CMs by Day')

plt.tight_layout()
plt.show()