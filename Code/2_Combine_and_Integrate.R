# Refactored script for combining and integrating the three cardiomyocyte datasets
# This script should be run from the project root directory.

library(Seurat)
library(dplyr)
library(ggplot2)
source("Code/utils.R") # Load the shared functions

# 1. Configuration
PATH_ADULT_CM <- "ProcessedData/AdultPrimary_Cardiomyocytes.rds"
PATH_FETAL_CM <- "ProcessedData/FetalPrimary_Cardiomyocytes.rds"
PATH_HIPSC_CM <- "ProcessedData/hipsc-cm_Cardiomyocytes.rds"
OUTPUT_INTEGRATED_RDS <- "ProcessedData/Integrated_Cardiomyocytes_All_Origins.rds"

# 2. Load Processed Seurat Objects
cat("Loading individual cardiomyocyte datasets...\n")
adult_cm <- readRDS(PATH_ADULT_CM)
fetal_cm <- readRDS(PATH_FETAL_CM)
hipsc_cm <- readRDS(PATH_HIPSC_CM)

# Clean up objects for merging
adult_cm <- DietSeurat(adult_cm, assays = "RNA", layers = "counts")
fetal_cm <- DietSeurat(fetal_cm, assays = "RNA", layers = "counts")
hipsc_cm <- DietSeurat(hipsc_cm, assays = "RNA", layers = "counts")
gc()

# 3. Integrate using the new function
# The function handles merging, normalization, scaling, PCA, and Harmony
all_cm_list <- list('Adult Primary' = adult_cm, 'Fetal Primary' = fetal_cm, 'hiPSC-CM' = hipsc_cm)

# The 'origin' column must exist, so let's add it from the names of the list
for (origin_name in names(all_cm_list)) {
  all_cm_list[[origin_name]]$origin <- origin_name
}

# The IntegrateSeuratObjects function expects a list of objects to merge
integrated_data <- IntegrateSeuratObjects(all_cm_list, group_by_vars = "origin", n_features = 3000, npcs = 40)

#if already saved, run this to skip all above, otherwise DO NOT RUN
integrated_data<- readRDS("ProcessedData/Integrated_Cardiomyocytes_All_Origins.rds")

# 4. Visualization and Saving
cat("Generating visualization plots...\n")
p1 <- DimPlot(integrated_data, group.by = "origin", pt.size = 0.5) +
  ggtitle("Integrated UMAP by Origin")
ggsave("figures/Seurat/Integrated_UMAP_by_Origin.png", plot = p1, width = 8, height = 6)

p2 <- DimPlot(integrated_data, group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  ggtitle("Integrated UMAP by Cluster")
ggsave("figures/Seurat/Integrated_UMAP_by_Cluster.png", plot = p2, width = 8, height = 6)

cat("Saving the final integrated Seurat object to", OUTPUT_INTEGRATED_RDS, "...\n")
saveRDS(integrated_data, OUTPUT_INTEGRATED_RDS)

p3 <- DimPlot(integrated_data, group.by = "cell_type")
ggsave("figures/Seurat/Integrated_UMAP_by_Cell_type.png", plot = p3, width = 8, height = 6)



cat("Analysis complete!\n")
