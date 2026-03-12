# Reusable functions for the BioMarker Identification project

#' Process a Seurat Object through a Standard Workflow
#'
#' This function takes a Seurat object and applies a standard single-sample
#' processing workflow including normalization, feature selection, scaling, and
#' dimensionality reduction.
#'
#' @param seu A Seurat object.
#' @param n_features The number of variable features to identify.
#' @param regress_vars A character vector of metadata variables to regress out during scaling.
#' @param npcs The number of principal components to compute.
#' @return A processed Seurat object.
#'
ProcessSeuratObject <- function(seu, n_features = 2000, regress_vars = NULL, npcs = 30) {
  
  cat("Normalizing data...\n")
  seu <- NormalizeData(seu, verbose = FALSE)
  
  cat("Finding variable features...\n")
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
  
  cat("Scaling data on variable features only...\n")
  seu <- ScaleData(seu, vars.to.regress = regress_vars, features = VariableFeatures(object = seu), verbose = FALSE)
  
  cat("Running PCA...\n")
  seu <- RunPCA(seu, npcs = npcs, verbose = FALSE)
  
  return(seu)
}


#' Integrate Seurat Objects using Harmony
#'
#' This function takes a list of Seurat objects or a merged Seurat object and
#' runs Harmony integration to correct for specified batch effects.
#'
#' @param data_object A merged Seurat object or a list of Seurat objects to be merged.
#' @param group_by_vars The metadata variable(s) to use for batch correction.
#' @param n_features The number of variable features to use for integration.
#' @param npcs The number of principal components to use for Harmony.
#' @return An integrated Seurat object.
#'
IntegrateSeuratObjects <- function(data_object, group_by_vars, n_features = 3000, npcs = 50) {
  
  # If the input is a list, merge it first
  if (is.list(data_object)) {
    cat("Merging", length(data_object), "Seurat objects...
")
    data_object <- merge(data_object[[1]], data_object[-1])
  }
  
  cat("Running standard processing before integration...
")
  data_object <- NormalizeData(data_object, verbose = FALSE)
  data_object <- FindVariableFeatures(data_object, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
  data_object <- ScaleData(data_object, features = VariableFeatures(object = data_object), verbose = FALSE)
  data_object <- RunPCA(data_object, npcs = npcs, verbose = FALSE)
  
  cat("Running Harmony integration on '", group_by_vars, "'...
")
  data_object <- RunHarmony(
    object = data_object,
    group.by.vars = group_by_vars,
    reduction = "pca",
    assay.use = "RNA",
    project.dim = FALSE
  )
  
  cat("Running UMAP and finding clusters on integrated data...
")
  data_object <- FindNeighbors(data_object, reduction = "harmony", dims = 1:npcs)
  data_object <- FindClusters(data_object, resolution = 0.5)
  data_object <- RunUMAP(data_object, reduction = "harmony", dims = 1:npcs)

  # In Seurat v5, it's good practice to JoinLayers after integration
  if ("harmony" %in% Reductions(data_object)) {
    data_object <- JoinLayers(data_object)
  }
  
  return(data_object)
}
