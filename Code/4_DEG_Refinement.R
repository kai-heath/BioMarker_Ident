library(tidyr)
library(dplyr)
library(biomaRt)


adultVFetal <- tibble(read.csv("Analysis_Results/adult_vs_fetal_primary_markers.csv"))
adultVHipsc <- tibble(read.csv("Analysis_Results/adult_vs_hipsc_markers.csv"))
FetalVHipsc <- tibble(read.csv("Analysis_Results/fetal_vs_hipsc_markers.csv"))

# 1. Merge the two lists on Gene Symbol
combined_metrics <- full_join(adultVFetal, adultVHipsc, by = "X", suffix = c("_Native", "_IPSC"))

# 2. Calculate the Power Score
combined_metrics$Maturation_Power <- combined_metrics$avg_log2FC_Native * combined_metrics$avg_log2FC_IPSC

# 3. Sort by the Score
combined_metrics <- combined_metrics[order(-combined_metrics$Maturation_Power), ]

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get biotypes for your 278 genes
gene_info <- getBM(attributes = c('external_gene_name', 'gene_biotype'),
                   filters = 'external_gene_name',
                   values = combined_metrics$X,
                   mart = ensembl)

coding_genes <- gene_info$external_gene_name[gene_info$gene_biotype == "protein_coding"]



final_protein_coding_list <- combined_metrics[combined_metrics$X %in% coding_genes, ] %>%
  filter(pct.1_Native > 0.0, pct.2_Native > 0.0, pct.1_IPSC > 0.0, pct.2_IPSC > 0.0)

final_protein_coding_list <- filter(final_protein_coding_list, p_val_adj_Native < 0.05)
final_protein_coding_list <- filter(final_protein_coding_list, p_val_adj_IPSC < 0.05)

write.csv(final_protein_coding_list, "Analysis_Results/DEGs_Merged.csv")

hipsc_cm <- readRDS("")