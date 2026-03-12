# Step 3: Marker Analysis
# This script is for downstream analysis of the integrated cardiomyocyte object.
# This script should be run from the project root directory.

library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Configuration
PATH_INTEGRATED_RDS <- "ProcessedData/Integrated_Cardiomyocytes_All_Origins.rds"
OUTPUT_DIR <- "Analysis_Results/"

# Create the output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

# 2. Load the Integrated Seurat Object
cat("Loading the integrated cardiomyocyte object...\n")
integrated_cm <- readRDS(PATH_INTEGRATED_RDS)

# 3. Perform Differential Gene Expression (DGE) Analysis
# This is where you will compare the different origins to find maturation markers.

cat("Finding markers between Adult Primary and Fetal Primary...\n")
Idents(integrated_cm) <- integrated_cm$origin
adult_vs_fetal_markers <- FindMarkers(integrated_cm, ident.1 = "Adult Primary", ident.2 = "Fetal Primary", logfc.threshold = 0.25, min.pct = 0.1)
write.csv(adult_vs_fetal_markers, file.path(OUTPUT_DIR, "adult_vs_fetal_primary_markers.csv"))

cat("Finding markers between Adult Primary and hiPSC-CM...\n")
adult_vs_hipsc_markers <- FindMarkers(integrated_cm, ident.1 = "Adult Primary", ident.2 = "hiPSC-CM", logfc.threshold = 0.25, min.pct = 0.1)
write.csv(adult_vs_hipsc_markers, file.path(OUTPUT_DIR, "adult_vs_hipsc_markers.csv"))

cat("Finding markers between Fetal Primary and hiPSC-CM...\n")
fetal_vs_hipsc_markers <- FindMarkers(integrated_cm, ident.1 = "Fetal Primary", ident.2 = "hiPSC-CM", logfc.threshold = 0.25, min.pct = 0.1)
write.csv(fetal_vs_hipsc_markers, file.path(OUTPUT_DIR, "fetal_vs_hipsc_markers.csv"))


# 4. Visualization
# Example: Volcano plot for Adult vs Fetal Primary
volcano_plot <- ggplot(adult_vs_fetal_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 1, "Significant", "Not Significant")), alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Adult vs. Fetal Primary Cardiomyocytes",
       x = "Average Log2 Fold Change",
       y = "-log10(Adjusted p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

ggsave(file.path(OUTPUT_DIR, "volcano_adult_vs_fetal_primary.png"), plot = volcano_plot, width = 8, height = 6)


cat("Marker analysis script is set up. You can add more complex analyses here.\n")
