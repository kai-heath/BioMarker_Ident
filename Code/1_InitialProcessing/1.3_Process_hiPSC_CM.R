library(Seurat)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)
source("Code/utils.R")

# grab GSE202398 file path
PATH_GSE202398 <- "RawData/GSE202398_(D0-D30Diff)/"

# Find only the 'filtered' h5 files
h5_files <- list.files(PATH_GSE202398, pattern = "_filtered_feature_bc_matrix\\.h5$", full.names = TRUE)

# this can be done multithreaded but I kept having issues so I swtiched to sequential
gse202398_list <- lapply(h5_files, function(h5_file) {

  sample_id <- gsub("_filtered_feature_bc_matrix\\.h5$", "", basename(h5_file))
  counts <- Read10X_h5(h5_file)
  
  # Filter HTOs that have very few total counts to avoid HTODemux error
  hto_counts <- counts$`Antibody Capture`
  valid_htos <- rownames(hto_counts)[rowSums(as.matrix(hto_counts)) > 1000]
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

  # HTODemux takes the tagged data and demultiplexes the data to identify the origin of each cell. In this case, what age/day the cell is and cell type
  seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR")
  seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.99)

  seu$sampleID <- sample_id
  seu$dataset <- "GSE202398" # Add a dataset identifier
  seu <- JoinLayers(seu)
  # Basic QC
  seu[["percent_mito"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent_mito < 20)

  seu <- NormalizeData(seu)
  # Find and scale variable features
  seu <- FindVariableFeatures(seu, selection.method = "mean.var.plot")
  seu <- ScaleData(seu, features = VariableFeatures(seu))
  
  # use HTO demux to remove doublets
  cat("Filtering for singlets based on HTO...\n")
  seu <- subset(seu, subset = HTO_classification.global == "Singlet")
  
  # Map HTO names to timepoint and cell_line
  seu$hash.ID <- as.character(seu$hash.ID)
  
  # Split hash.ID to get cell_line and day
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

# Combine the two datasets and Integrate
# To merge, they must have the same assays and layers.
gse202398_merged <- DietSeurat(gse202398_merged, assays = "RNA", layers = "counts")
gse263372_seu <- DietSeurat(gse263372_seu, assays = "RNA", layers = "counts")
gc()

# Use utils function to integrate
gc()
hipsc_cm_integrated <- IntegrateSeuratObjects(
  c(gse202398_merged, gse263372_seu),
  group_by_vars = "dataset", # Correct for the two different sources
  n_features = 3000,
  npcs = 40
)

# Identify Cardiomyocytes
cardiomyocyte_markers <- list(c("TNNT2", "MYH7", "MYL2", "MYL7", "ACTN2"))
hipsc_cm_integrated <- AddModuleScore(hipsc_cm_integrated, features = cardiomyocyte_markers, name = "CM_Score")

# Use a quantile threshold to identify CMs
cm_threshold <- quantile(hipsc_cm_integrated$CM_Score1, 0.30)
cardiomyocytes <- subset(hipsc_cm_integrated, subset = CM_Score1 > cm_threshold)

# add  metadata
hipsc_cm_integrated$origin <- "hiPSC-CM"
cardiomyocytes$origin <- "hiPSC-CM"

# Define output paths
OUTPUT_CM_RDS <- "ProcessedData/hipsc-cm_Cardiomyocytes.rds"
OUTPUT_FULL_RDS <- "ProcessedData/hipsc-cm_Integrated.rds"

# Save both the final integrated object and the subsetted CM object
saveRDS(cardiomyocytes, OUTPUT_CM_RDS)
saveRDS(hipsc_cm_integrated, OUTPUT_FULL_RDS)

# Generate diagnostic plots if wanted/needed
p1 <- DimPlot(hipsc_cm_integrated, group.by = "dataset", pt.size = 0.5) + ggtitle("hiPSC-CM by Dataset")
p2 <- DimPlot(hipsc_cm_integrated, group.by = "timepoint", pt.size = 0.5, label = TRUE) + ggtitle("hiPSC-CM by Timepoint")
p3 <- FeaturePlot(hipsc_cm_integrated, "CM_Score1") + ggtitle("Cardiomyocyte Score")

# Save plots
ggsave("Integrated_UMAP_by_Dataset_hiPSC_CM.png", plot = p1, width = 8, height = 6)
ggsave("Integrated_UMAP_by_Timepoint_hiPSC_CM.png", plot = p2, width = 10, height = 6)
ggsave("Integrated_UMAP_CM_Score_hiPSC_CM.png", plot = p3, width = 8, height = 6)

p2