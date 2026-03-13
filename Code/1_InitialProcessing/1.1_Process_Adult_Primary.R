# Refactored script for processing Adult Primary data using shared utility functions
# This script should be run from the project root directory.

library(Seurat)
library(Matrix)
library(rhdf5)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
source("Code/utils.R") # Load the shared functions

# 1. Configuration
PATH_H5AD <- "RawData/HCA_Adult_Primary_Tissue/Global_raw.h5ad"
OUTPUT_RDS <- "ProcessedData/AdultPrimary_Cardiomyocytes.rds"

# 2. Surgical Extraction of Cardiomyocytes (this part is unique and remains)
cat("Reading metadata to find cardiomyocyte indices...\n")
obs_cell_types <- as.numeric(h5read(PATH_H5AD, "/obs/cell_type"))
barcodes <- as.character(h5read(PATH_H5AD, "/obs/barcode"))
cm_indices <- which(obs_cell_types %in% c(0, 1)) # Atrial and Ventricular CMs
cat("Found", length(cm_indices), "cardiomyocytes.\n")

cat("Extracting sparse matrix data for CMs only...\n")
indptr <- as.numeric(h5read(PATH_H5AD, "/X/indptr"))
genes <- as.character(h5read(PATH_H5AD, "/var/_index"))
all_indices <- h5read(PATH_H5AD, "/X/indices")
all_data <- h5read(PATH_H5AD, "/X/data")

start_indices <- indptr[cm_indices] + 1
end_indices <- indptr[cm_indices + 1]
extract_idx <- unlist(mapply(seq, start_indices, end_indices, SIMPLIFY = FALSE))

subset_indices <- all_indices[extract_idx]
subset_data <- as.numeric(all_data[extract_idx])
rm(all_indices, all_data); gc()

row_lengths <- end_indices - start_indices + 1
new_indptr <- c(0, cumsum(row_lengths))

counts_csr <- sparseMatrix(
  j = subset_indices + 1,
  p = new_indptr,
  x = subset_data,
  dims = c(length(cm_indices), length(genes))
)
rownames(counts_csr) <- barcodes[cm_indices]
colnames(counts_csr) <- genes
counts <- t(counts_csr)
rm(counts_csr); gc()

# 3. Create Seurat Object and Add Metadata
cat("Creating Seurat object...\n")
seu <- CreateSeuratObject(counts = counts)
rm(counts); gc()

# Add metadata from the H5AD file
cell_labels <- as.character(h5read(PATH_H5AD, "/obs/__categories/cell_type"))
donor_ids <- as.character(h5read(PATH_H5AD, "/obs/__categories/donor"))[as.numeric(h5read(PATH_H5AD, "/obs/donor")) + 1]
region <- as.character(h5read(PATH_H5AD, "/obs/__categories/region"))[as.numeric(h5read(PATH_H5AD, "/obs/region")) + 1]

seu$cell_type <- cell_labels[obs_cell_types[cm_indices] + 1]
seu$donor <- donor_ids[cm_indices]
seu$region <- region[cm_indices]
seu$percent_mito <- as.numeric(h5read(PATH_H5AD, "/obs/pct_counts_mt"))[cm_indices]
seu$origin <- "Adult Primary"
gc()

# 4. Doublet Detection
cat("Running scDblFinder...\n")
sce <- as.SingleCellExperiment(seu)
sce <- scDblFinder(sce)
seu$scDblFinder.class <- sce$scDblFinder.class
seu <- subset(seu, subset = scDblFinder.class != "doublet")
cat("Removed potential doublets.\n")

# 5. Process using the new reusable function
# This replaces NormalizeData, FindVariableFeatures, ScaleData, RunPCA
rm(counts, extract_idx, start_indices, end_indices, row_lengths, new_indptr, subset_indices, subset_data, barcodes, genes, obs_cell_types, cm_indices, donor_ids, region, cell_labels); gc()
cat("Processing Seurat object with standard workflow...\n")
seu <- ProcessSeuratObject(seu, n_features = 2000, npcs = 30)

# The function doesn't do clustering/UMAP, so we do it here if needed
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)

# 6. Save
cat("Saving to", OUTPUT_RDS, "...\n")
saveRDS(seu, OUTPUT_RDS)
cat("Refactored Adult Primary processing complete!\n")
