# Script to process the two hiPSC-CM datasets (GSE202398 and GSE263372)
# This script should be run from the project root directory.

library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)
library(scDblFinder)
library(SingleCellExperiment)
source("Code/utils.R")

# ---
# Part 1: Process GSE202398 (D0-D30 Diff) from H5 files
# ---
cat("--- Part 1: Processing GSE202398 ---\n")
PATH_GSE202398 <- "RawData/GSE202398_(D0-D30Diff)/"

# Set parallel plan for loading
plan("multisession", workers = 4)

# Find only the 'filtered' h5 files
h5_files <- list.files(PATH_GSE202398, pattern = "_filtered_feature_bc_matrix\\.h5$", full.names = TRUE)

gse202398_list <- lapply(h5_files, function(h5_file) {
  sample_id <- gsub("_filtered_feature_bc_matrix\\.h5$", "", basename(h5_file))
  cat("Loading sample:", sample_id, "\n")
  
  counts <- Read10X_h5(h5_file)
  seu <- CreateSeuratObject(counts = counts)
  
  seu$sampleID <- sample_id
  seu$dataset <- "GSE202398" # Add a dataset identifier
  seu <- JoinLayers(seu)
  # Basic QC
  seu[["percent_mito"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent_mito < 20)
  
  # Doublet detection
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce)
  seu$scDblFinder.class <- sce$scDblFinder.class
  seu <- subset(seu, subset = scDblFinder.class != "doublet")
  
  if (ncol(seu) < 100) return(NULL)
  return(seu)
})

plan(sequential) # Return to sequential processing

# Filter out any NULLs and merge
gse202398_list <- Filter(Negate(is.null), gse202398_list)
if (length(gse202398_list) == 0) {
    stop("No valid samples were processed from GSE202398.")
}

cat("Merging", length(gse202398_list), "samples from GSE202398...\n")
if (length(gse202398_list) > 1) {
  gse202398_merged <- merge(gse202398_list[[1]], gse202398_list[-1])
} else {
  gse202398_merged <- gse202398_list[[1]]
}
rm(gse202398_list); gc()

# ---
# Part 2: Load GSE263372 (BioReactor) from RDS file
# ---
cat("\n--- Part 2: Loading GSE263372 ---\n")
PATH_GSE263372_RDS <- "RawData/GSE263372_(BioReator)/GSE263372_combined_seurat.rds"

gse263372_seu <- readRDS(PATH_GSE263372_RDS)
cat("Loaded GSE263372 object with", ncol(gse263372_seu), "cells.\n")

# Add a dataset identifier. We need to check existing metadata to avoid conflicts.
gse263372_seu$dataset <- "GSE263372"

# Perform a quick QC check on the loaded object
gse263372_seu[["percent_mito"]] <- PercentageFeatureSet(gse263372_seu, pattern = "^MT-")
gse263372_seu <- subset(gse263372_seu, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent_mito < 20)
sce <- as.SingleCellExperiment(gse263372_seu)
sce <- scDblFinder(sce)
gse263372_seu$scDblFinder.class <- sce$scDblFinder.class
gse263372_seu <- subset(gse263372_seu, subset = scDblFinder.class != "doublet")
cat("After QC and doublet removal, GSE263372 object has", ncol(gse263372_seu), "cells.\n")


# ---
# Part 3: Combine the two datasets and Integrate
# ---
cat("\n--- Part 3: Combining and Integrating Datasets ---\n")

# To merge, they must have the same assays and layers.
# Let's simplify them to just the 'RNA' assay with 'counts'.
gse202398_merged <- DietSeurat(gse202398_merged, assays = "RNA", layers = "counts")
gse263372_seu <- DietSeurat(gse263372_seu, assays = "RNA", layers = "counts")
gc()


# Use our utility function to integrate them, correcting for the 'dataset' effect
gc()
hipsc_cm_integrated <- IntegrateSeuratObjects(
  c(gse202398_merged, gse263372_seu),
  group_by_vars = "dataset", # Correct for the two different sources
  n_features = 3000,
  npcs = 40
)

# ---
# Part 4: Identify Cardiomyocytes and Save
# ---
cat("\n--- Part 4: Identifying Cardiomyocytes ---\n")
cardiomyocyte_markers <- list(c("TNNT2", "MYH7", "MYL2", "MYL7", "ACTN2"))
hipsc_cm_integrated <- AddModuleScore(hipsc_cm_integrated, features = cardiomyocyte_markers, name = "CM_Score")

# Use a quantile threshold to identify CMs
cm_threshold <- quantile(hipsc_cm_integrated$CM_Score1, 0.30) # Use a slightly more stringent threshold
cat("Using CM Score threshold:", cm_threshold, "\n")
cardiomyocytes <- subset(hipsc_cm_integrated, subset = CM_Score1 > cm_threshold)
cat("Found", ncol(cardiomyocytes), "cardiomyocytes out of", ncol(hipsc_cm_integrated), "total cells.\n")

# Final step: add the 'origin' metadata for the final merge
hipsc_cm_integrated$origin <- "hiPSC-CM"
cardiomyocytes$origin <- "hiPSC-CM"

# Define output paths
OUTPUT_CM_RDS <- "ProcessedData/hipsc-cm_Cardiomyocytes.rds"
OUTPUT_FULL_RDS <- "ProcessedData/hipsc-cm_Integrated.rds"

# Save both the final integrated object and the subsetted CM object
cat("Saving final hiPSC-CM cardiomyocyte object to", OUTPUT_CM_RDS, "\n")
saveRDS(cardiomyocytes, OUTPUT_CM_RDS)

cat("Saving final hiPSC-CM integrated (all cells) object to", OUTPUT_FULL_RDS, "\n")
saveRDS(hipsc_cm_integrated, OUTPUT_FULL_RDS)

cat("hiPSC-CM processing complete!\n")

# Generate a diagnostic plot
DimPlot(hipsc_cm_integrated, group.by = "stim", pt.size = 0.5)

FeaturePlot(hipsc_cm_integrated, "CM_Score1")

ggsave("hipsc_cm_integrated_by_dataset.png", plot = p1, width = 8, height = 6)
rm(p1)

