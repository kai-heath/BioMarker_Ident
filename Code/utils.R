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

#' Integrate Seurat Objects using a Memory-Efficient Harmony Workflow
#'
#' This function takes a list of Seurat objects and integrates them using a
#' memory-optimized workflow. It processes each object individually before merging,
#' selects common features, and scales only those features on the final merged object.
#'
#' @param object_list A list of Seurat objects to integrate.
#' @param group_by_vars The metadata variable(s) to use for batch correction.
#' @param n_features The number of variable features to identify and use for integration.
#' @param npcs The number of principal components to compute and use for Harmony.
#' @return A single, integrated Seurat object.
IntegrateSeuratObjects <- function(object_list, group_by_vars, n_features = 2000, npcs = 30) {
  library(harmony)

  # Handle the case of a single object gracefully
  if (length(object_list) == 1) {
    cat("Only one object found, running standard processing instead of integration...\n")
    data_object <- ProcessSeuratObject(object_list[[1]], n_features = n_features, npcs = npcs)
    # Downstream analysis for single object
    data_object <- FindNeighbors(data_object, dims = 1:npcs)
    data_object <- FindClusters(data_object, resolution = 0.5)
    data_object <- RunUMAP(data_object, dims = 1:npcs)
    return(data_object)
  }

  cat("Step 1: Normalizing and finding variable features for each of the", length(object_list), "objects...\n")
  object_list <- lapply(object_list, function(obj) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
    return(obj)
  })

  cat("Step 2: Selecting common features for integration...\n")
  integration_features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = n_features)
  
  cat("Step 3: Merging objects...\n")
  merged_object <- merge(object_list[[1]], y = object_list[-1])
  rm(object_list); gc()

  cat("Step 4: Scaling and PCA on MERGED data using common features only...\n")
  # This is the key memory-saving step.
  merged_object <- ScaleData(merged_object, features = integration_features, verbose = FALSE)
  merged_object <- RunPCA(merged_object, features = integration_features, npcs = npcs, verbose = FALSE)
  
  cat("Step 5: Running Harmony integration...\n")
  merged_object <- RunHarmony(
    object = merged_object,
    group.by.vars = group_by_vars,
    assay.use = "RNA" # Explicitly use the RNA assay
  )
  
  cat("Step 6: Downstream analysis (UMAP, Clustering)...\n")
  merged_object <- FindNeighbors(merged_object, reduction = "harmony", dims = 1:npcs)
  merged_object <- FindClusters(merged_object, resolution = 0.5)
  merged_object <- RunUMAP(merged_object, reduction = "harmony", dims = 1:npcs)

  # Seurat v5 post-integration step
  # if ("harmony" %in% Reductions(merged_object)) {
  #   merged_object <- JoinLayers(merged_object)
  # }
  
  return(merged_object)
}
