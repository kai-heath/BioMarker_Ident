library(Seurat)

# runs a bunch of functions needed for each seurat object made.
#' @param seu A Seurat object.
#' @param n_features The number of variable features to identify.
#' @param regress_vars A character vector of metadata variables to regress out during scaling.
#' @param npcs The number of principal components to compute.
#' @return A processed Seurat object.
#'
ProcessSeuratObject <- function(seu, n_features = 2000, regress_vars = NULL, npcs = 30) {
  
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
  seu <- ScaleData(seu, vars.to.regress = regress_vars, features = VariableFeatures(object = seu), verbose = FALSE)
  seu <- RunPCA(seu, npcs = npcs, verbose = FALSE)
  
  return(seu)
}


# combines seurat objects using harmony
#' @param object_list A list of Seurat objects to integrate.
#' @param group_by_vars The metadata variable(s) to use for batch correction.
#' @param n_features The number of variable features to identify and use for integration.
#' @param npcs The number of principal components to compute and use for Harmony.
#' @return A single, integrated Seurat object.
IntegrateSeuratObjects <- function(object_list, group_by_vars, n_features = 2000, npcs = 30) {
  library(harmony)

  # If for whatever reason theres just one object, this catches that
  if (length(object_list) == 1) {
    data_object <- ProcessSeuratObject(object_list[[1]], n_features = n_features, npcs = npcs)
    data_object <- FindNeighbors(data_object, dims = 1:npcs)
    data_object <- FindClusters(data_object, resolution = 0.5)
    data_object <- RunUMAP(data_object, dims = 1:npcs)
    return(data_object)
  }

  # Ideally the above function was ran on the data beforehand, so they'll have all that done already, but just in case, some critical parts are rerun (doesn't cause issues)
  object_list <- lapply(object_list, function(obj) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
    return(obj)
  })

  integration_features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = n_features)
  

  merged_object <- merge(object_list[[1]], y = object_list[-1])
  rm(object_list); gc()


  merged_object <- ScaleData(merged_object, features = integration_features, verbose = FALSE)
  merged_object <- RunPCA(merged_object, features = integration_features, npcs = npcs, verbose = FALSE)
  
  merged_object <- RunHarmony(
    object = merged_object,
    group.by.vars = group_by_vars,
    assay.use = "RNA" # Explicitly use the RNA assay
  )

  merged_object <- FindNeighbors(merged_object, reduction = "harmony", dims = 1:npcs)
  merged_object <- FindClusters(merged_object, resolution = 0.5)
  merged_object <- RunUMAP(merged_object, reduction = "harmony", dims = 1:npcs)


   if ("harmony" %in% Reductions(merged_object)) {
     merged_object <- JoinLayers(merged_object)
   }
  
  return(merged_object)
}
