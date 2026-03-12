library(Seurat)
library(Matrix)
library(future)
library(RcppHNSW)
library(igraph)
library(future.apply)
library(niceRplots)
library("scDblFinder")
library(ggplot2)
library(ROCR)
library(KernSmooth)
library(dplyr)
library(patchwork)
library(harmony)
library(gprofiler2)
library(plotly)

#setting up path to data vars
PATH_DATA <- "RawData/FetalTissue/1_Chromium_cellranger_data_SC/"
file_list <- list.files(PATH_DATA, full.names = FALSE, pattern = ".h5")

# load the matrix, create a Seurat object
# and store it in a list for all datasets in parallel
data_list <- future_lapply(file_list, function(x) {
  
  # 1. Load Matrix & Initialize Object
  counts <- Matrix::Matrix(rhdf5::h5read(paste0(PATH_DATA, x), name = "/shoji/Expression"), sparse = TRUE)
  cell_ids <- as.character(unlist(rhdf5::h5read(paste0(PATH_DATA, x), name = "/shoji/Cellid")))
  gene_names <- as.character(unlist(rhdf5::h5read(paste0(PATH_DATA, x), name = "/shoji/Gene")))

  rownames(counts) <- make.unique(as.character(gene_names))
  colnames(counts) <- as.character(cell_ids)
  
  seu <- CreateSeuratObject(counts = counts)
  seu$sampleID <- gsub("\\.h5$", "", x)
  
  # 2. Add Pre-computed H5 Metadata 
  # Using data.frame to ensure clean rowname matching with Seurat AddMetaData
  meta_df <- data.frame(
    Shoji_DoubletFlag = as.character(unlist(rhdf5::h5read(paste0(PATH_DATA, x), name = "/shoji/DoubletFlag"))),
    Shoji_DoubletScore = as.numeric(unlist(rhdf5::h5read(paste0(PATH_DATA, x), name = "/shoji/DoubletScore"))),
    Shoji_DropletClass = as.character(unlist(rhdf5::h5read(paste0(PATH_DATA, x), name = "/shoji/DropletClass"))),
    row.names = cell_ids
  )
  seu <- AddMetaData(seu, metadata = meta_df)
  
  # 3. Calculate QC Metrics
  seu <- PercentageFeatureSet(seu, pattern = "^RP[LS]", col.name = "percent_ribo")
  seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent_mito")
  seu <- PercentageFeatureSet(seu, pattern = "^HB[^ESP]|HBE1", col.name = "percent_hb")
  seu <- PercentageFeatureSet(seu, pattern = "^HSP", col.name = "percent_hsp")
  seu$apoptotic_maybe <- seu$percent_hsp > 2
  
  # 4. Filter Cells (Using your established logic)
  cells_to_filter <- (seu$percent_mito > 30) | 
                     (seu$percent_ribo < 3 & seu$percent_hb < 10) | 
                     (seu$nFeature_RNA < 250 & seu$percent_hb < 10)
  
  seu <- seu[, !cells_to_filter]
  
  # Failsafe: If a sample is mostly dead cells and filters down to <50 cells, 
  # PCA/DoubletFinder will crash the whole loop. Skip it safely.
  if(ncol(seu) < 50) return(NULL)
  
  # 5. Preprocess specifically for DoubletFinder
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
  
  sce <- as.SingleCellExperiment(seu)
  
  # F. Run scDblFinder
  # It automatically normalizes, runs PCA, calculates expected rates, and scores doublets
  sce <- scDblFinder(sce)
  
  # G. Pull the classifications back into your Seurat object
  seu$DF_classification <- sce$scDblFinder.class # Returns "singlet" or "doublet"
  seu$doublet_score <- sce$scDblFinder.score
  
  # 7. Memory Cleanup 
  seu <- DietSeurat(seu, layers = c("counts", "data", "scale.data"), dimreducs = NULL)
  
  return(seu)
})

# Clean up any NULL entries (samples that failed QC completely)
data_list <- Filter(Negate(is.null), data_list)
names(data_list) <- sapply(data_list, function(x) unique(x$sampleID))

# Force garbage collection
gc()


data <- merge(data_list[[1]], data_list[-1])
rm(data_list)
data <- subset(data, subset = DF_classification == "singlet")
gc()

plan(sequential)
data <- NormalizeData(data)
data <- CellCycleScoring(data, 
                          s.features = cc.genes.updated.2019$s.genes,
                          g2m.features = cc.genes.updated.2019$g2m.genes)

gc()

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data, npcs = 30, verbose = FALSE)

# 3. Run Harmony Integration using v5 IntegrateLayers
# This looks at the different sample layers and aligns them in PCA space
data <- IntegrateLayers(
  object = data, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = TRUE
)

# 4. Downstream Clustering and UMAP using the batch-corrected Harmony embeddings
data <- FindNeighbors(data, reduction = "harmony", dims = 1:30)
data <- FindClusters(data, resolution = 0.5)
data <- RunUMAP(data, reduction = "harmony", dims = 1:30)

DimPlot(data, group.by = "seurat_clusters", label = TRUE, reduction = "umap") + NoLegend()

data <- JoinLayers(data)

cardiomyocyte_markers <- c("TNNT2", "MYH7", "MYL2", "MYL7")

# Plot them on the UMAP
FeaturePlot(data, 
            features = cardiomyocyte_markers, 
            ncol = 2, 
            pt.size = 0.5, 
            label = TRUE, # Keeps the cluster numbers visible
            repel = TRUE)

# Extract the cardiomyocyte lineage
cardiomyocytes <- subset(data, idents = c(0, 6, 10, 13, 26))

saveRDS(data, "ProcessedData/FetalTissue_SeuratObject.rds")

saveRDS(cardiomyocytes, "ProcessedData/FetalTissue_Cardiomyocytes.rds")

DimPlot(cardiomyocytes, group.by = "seurat_clusters", label = TRUE, reduction = "umap") + NoLegend()
