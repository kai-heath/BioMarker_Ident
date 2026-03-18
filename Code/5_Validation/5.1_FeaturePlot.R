library(Seurat)
library(dplyr)
library(ggplot2)

integrated_cm <- readRDS("ProcessedData/hipsc-cm_Cardiomyocytes.rds")

# Define the gene sets based on your refined DEGs
# Using top consistent markers from your merged CSV
mature_genes <- c("MYOM3", "AQP7", "TMEM178B", "PXDNL", "CFAP61", 
                       "ASB18", "TOGARAM2", "UGT2B4", "RPL3L", "NRAP")

immature_genes <- list(c("MDK", "PTTG1", "KRT18", "MARCKSL1", "TNNI1", 
                         "KRT19", "VIM", "COL3A1", "STMN1", "COL1A2"))
                         
# 1. Calculate Mature Module Score
integrated_cm <- AddModuleScore(
  object = integrated_cm,
  features = mature_genes,
  name = "Mature_Score"
)

# 2. Calculate Immature Module Score
integrated_cm <- AddModuleScore(
  object = integrated_cm,
  features = immature_genes,
  name = "Immature_Score"
)

# 3. Create the final Maturation Index (Mature - Immature)
# Seurat appends a '1' to the name provided in AddModuleScore
integrated_cm$Maturation_Index <- integrated_cm$Mature_Score1 - integrated_cm$Immature_Score1

FeaturePlot(
  integrated_cm, 
  features = "Maturation_Index", 
  reduction = "umap",
  ncol = 2,
  order = TRUE,       # Plots cells with higher expression on top so they aren't hidden
  min.cutoff = 'q9',  # Sets the minimum color to the 9th percentile to reduce background noise
  label = TRUE        # Adds cluster labels to help orient where the 'mature' cells are
)

FeaturePlot(
  integrated_cm, 
  features = "Mature_Score1", 
  reduction = "umap",
  ncol = 2,
  order = TRUE,       # Plots cells with higher expression on top so they aren't hidden
  min.cutoff = 'q9',  # Sets the minimum color to the 9th percentile to reduce background noise
  label = TRUE        # Adds cluster labels to help orient where the 'mature' cells are
)

# 1. Group the continuous index into 5 sequential "Phases" of maturation
integrated_cm$Maturation_Phase <- cut(
  integrated_cm$Maturation_Index, 
  breaks = 5, 
  labels = c("Phase 1 (Immature)", "Phase 2", "Phase 3", "Phase 4", "Phase 5 (Mature)")
)

# 2. Set the active identity to these new phases
Idents(integrated_cm) <- "Maturation_Phase"

# 3. Define the genes you want to validate (combine your top markers)
genes_to_validate <- c(
  "MDK", "PTTG1", "TNNI1", "VIM",       # Expected to fade out
  "CDIN1", "MYOM3", "AQP7", "PXDNL"      # Expected to fade in
)

integrated_cm <- ScaleData(integrated_cm, features = genes_to_validate)

# 4. Generate the Heatmap ordered by phase
DoHeatmap(integrated_cm, features = genes_to_validate, group.by = "Maturation_Phase") +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3")) +
  ggtitle("Gene Expression Flow Across Maturation Phases")

DimPlot(integrated_cm, group.by = "Day")
