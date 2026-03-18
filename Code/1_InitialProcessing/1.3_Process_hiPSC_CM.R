# Script to process the two hiPSC-CM datasets (GSE202398 and GSE263372)
# This script should be run from the project root directory.

library(Seurat)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)
source("Code/utils.R")

# ---
# Part 1: Process GSE202398 (D0-D30 Diff) from H5 files
# ---
cat("--- Part 1: Processing GSE202398 ---\n")
PATH_GSE202398 <- "RawData/GSE202398_(D0-D30Diff)/"

# Find only the 'filtered' h5 files
h5_files <- list.files(PATH_GSE202398, pattern = "_filtered_feature_bc_matrix\\.h5$", full.names = TRUE)

gse202398_list <- lapply(h5_files, function(h5_file) {
  sample_id <- gsub("_filtered_feature_bc_matrix\\.h5$", "", basename(h5_file))
  cat("Loading sample:", sample_id, "\n")
  
  counts <- Read10X_h5(h5_file)
  
  # Filter HTOs that have very few total counts to avoid HTODemux error
  # (some HTOs are present in the H5 header but were not used in that specific Run)
  hto_counts <- counts$`Antibody Capture`
  valid_htos <- rownames(hto_counts)[rowSums(as.matrix(hto_counts)) > 1000]
  cat("Valid HTOs for sample", sample_id, ":", paste(valid_htos, collapse = ", "), "\n")
  hto_counts <- hto_counts[valid_htos, , drop = FALSE]
  
  HTO_assay <- CreateAssay5Object(counts = hto_counts)
  seu <- CreateSeuratObject(counts = counts)
  seu[["HTO"]] <- HTO_assay

  # Filter out cells that have ZERO total HTO counts across all valid tags
  # HTODemux fails if there are cells with no HTO reads at all
  hto_totals_per_cell <- colSums(LayerData(seu, assay = "HTO", layer = "counts"))
  RNA_totals_per_cell <- colSums(LayerData(seu, assay = "RNA", layer = "counts"))

  print(table(hto_totals_per_cell == 0))
  seu <- subset(seu, cells = names(hto_totals_per_cell[hto_totals_per_cell > 0]))
  seu <- subset(seu, cells = names(RNA_totals_per_cell[RNA_totals_per_cell > 0]))
  print(table(colSums(LayerData(seu, assay = "HTO", layer = "counts")) == 0))

  seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR")
  seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.99)

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


  seu <- NormalizeData(seu)
  # Find and scale variable features
  seu <- FindVariableFeatures(seu, selection.method = "mean.var.plot")
  seu <- ScaleData(seu, features = VariableFeatures(seu))
  
  # After HTODemux, we want singlets
  cat("Filtering for singlets based on HTO...\n")
  seu <- subset(seu, subset = HTO_classification.global == "Singlet")
  
  # Map HTO names to timepoint and cell_line
  # hash.ID is like "LMNA-DAY1" after Seurat auto-correction of "_" to "-"
  seu$hash.ID <- as.character(seu$hash.ID)
  
  # Split the hash.ID to get cell_line and day
  parts <- strsplit(seu$hash.ID, "-")
  seu$cell_line <- sapply(parts, `[`, 1)
  seu$timepoint <- sapply(parts, `[`, 2)
  
  # Standardize timepoint: "DAY1" -> "D1"
  seu$timepoint <- gsub("DAY", "D", seu$timepoint)
  
  if (ncol(seu) < 100) return(NULL)
  return(seu)
})

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
# Map stim to timepoint for consistency when merging
gse263372_seu$timepoint <- gse263372_seu$stim

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

# Generate diagnostic plots
p1 <- DimPlot(hipsc_cm_integrated, group.by = "dataset", pt.size = 0.5) + ggtitle("hiPSC-CM by Dataset")
p2 <- DimPlot(hipsc_cm_integrated, group.by = "timepoint", pt.size = 0.5, label = TRUE) + ggtitle("hiPSC-CM by Timepoint")
p3 <- FeaturePlot(hipsc_cm_integrated, "CM_Score1") + ggtitle("Cardiomyocyte Score")

# Save plots
ggsave("Integrated_UMAP_by_Dataset_hiPSC_CM.png", plot = p1, width = 8, height = 6)
ggsave("Integrated_UMAP_by_Timepoint_hiPSC_CM.png", plot = p2, width = 10, height = 6)
ggsave("Integrated_UMAP_CM_Score_hiPSC_CM.png", plot = p3, width = 8, height = 6)

cat("hiPSC-CM processing complete!\n")
p2

