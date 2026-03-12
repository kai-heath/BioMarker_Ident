# Refactored script for processing Fetal Primary data using shared utility functions
# This script should be run from the project root directory.

library(Seurat)
library(Matrix)
library(future)
library(future.apply)
library(dplyr)
library(rhdf5)
source("Code/utils.R") # Load the shared functions

# 1. Configuration
PATH_DATA <- "RawData/FetalTissue/1_Chromium_cellranger_data_SC/"
OUTPUT_RDS <- "ProcessedData/FetalPrimary_SeuratObject.rds"
OUTPUT_CM_RDS <- "ProcessedData/FetalPrimary_Cardiomyocytes.rds"

# Set parallel plan
plan("multisession", workers = 4)

# 2. Sample Discovery and Loading
file_list <- list.files(PATH_DATA, full.names = FALSE, pattern = ".h5")

data_list <- future_lapply(file_list, function(x) {
    h5_file <- paste0(PATH_DATA, x)
    sample_id <- gsub("\.h5$", "", x)
    
    # Filter for Heart tissue early
    tissues <- as.character(unlist(rhdf5::h5read(h5_file, name = "/shoji/Tissue")))
    indices <- which(tissues == "Heart")
    if(length(indices) == 0) return(NULL)
    
    # Load Matrix and create Seurat object (same as original script)
    counts <- Matrix::Matrix(rhdf5::h5read(h5_file, name = "/shoji/Expression", index = list(NULL, indices)), sparse = TRUE)
    cell_ids <- as.character(unlist(rhdf5::h5read(h5_file, name = "/shoji/Cellid", index = list(indices))))
    gene_names <- as.character(unlist(rhdf5::h5read(h5_file, name = "/shoji/Gene")))
    rownames(counts) <- make.unique(as.character(gene_names))
    colnames(counts) <- as.character(cell_ids)
    
    seu <- CreateSeuratObject(counts = counts)
    seu$sampleID <- sample_id
    
    # QC and filtering (same as original script)
    seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_RNA > 250 & percent_mito < 30)
    
    if(ncol(seu) < 50) return(NULL)
    return(seu)
})

plan(sequential)
data_list <- Filter(Negate(is.null), data_list)
if(length(data_list) == 0) stop("No valid samples processed.")
names(data_list) <- sapply(data_list, function(x) unique(x$sampleID))
gc()

# 3. Process and Integrate using the new function
gc()
if (length(data_list) > 1) {
    data <- IntegrateSeuratObjects(data_list, group_by_vars = "sampleID", npcs = 30)
} else {
    data <- ProcessSeuratObject(data_list[[1]], npcs = 30)
    data <- FindNeighbors(data, dims = 1:30)
    data <- FindClusters(data, resolution = 0.5)
    data <- RunUMAP(data, dims = 1:30)
}


# 4. Robust CM Identification
cat("Identifying Cardiomyocytes...\n")
cardiomyocyte_markers <- list(c("TNNT2", "MYH7", "MYL2", "MYL7"))
data <- AddModuleScore(data, features = cardiomyocyte_markers, name = "CM_Score")

cm_check <- data@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(avg_score = mean(CM_Score1))
  
threshold <- median(cm_check$avg_score) + (sd(cm_check$avg_score) * 1.5)
cm_clusters <- cm_check %>% filter(avg_score > threshold) %>% pull(seurat_clusters)

cardiomyocytes <- subset(data, idents = cm_clusters)
cat("Found", ncol(cardiomyocytes), "cardiomyocytes.\n")

# 5. Add Origin Metadata and Save
data$origin <- "Fetal Primary"
cardiomyocytes$origin <- "Fetal Primary"

saveRDS(data, OUTPUT_RDS)
saveRDS(cardiomyocytes, OUTPUT_CM_RDS)
cat("Fetal Primary processing complete!\n")
