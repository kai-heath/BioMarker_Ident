library(Seurat)
library(Matrix)

# 1. Load your object
seurat_obj <- readRDS("RawData/GSE263372_(BioReator)/GSE263372_combined_seurat.rds")

# 2. Export the Count Matrix (Sparse Matrix format is best)
# We use .mtx because a standard .csv for 20k+ cells is too large for RAM
writeMM(obj = GetAssayData(seurat_obj, assay = "RNA", layer = "counts"), 
        file = "matrix.mtx")

# 3. Export Cell Metadata (obs)
write.csv(seurat_obj@meta.data, file = "metadata.csv", row.names = TRUE)

# 4. Export Gene Names (var)
# Ensure the order matches the matrix rows
write.csv(rownames(seurat_obj), file = "genes.csv", row.names = FALSE)
