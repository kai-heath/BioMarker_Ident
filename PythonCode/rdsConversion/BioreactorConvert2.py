import scanpy as sc
import pandas as pd

# 1. Load the matrix and transpose
adata = sc.read_mtx("RawData/GSE263372_(BioReator)/matrix.mtx").T

# 2. Load and assign gene names (var)
genes = pd.read_csv("RawData/GSE263372_(BioReator)/genes.csv")
adata.var_names = genes.iloc[:, 0].values

# 3. Load and assign cell metadata (obs)
metadata = pd.read_csv("RawData/GSE263372_(BioReator)/metadata.csv", index_col=0)
adata.obs = metadata

# 5. Final check and save
adata.write("RawData/GSE263372_(BioReator)/converted_data.h5ad")
