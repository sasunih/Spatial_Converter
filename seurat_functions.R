library("Seurat")
library("readr")
library("dplyr")
library("stringr")
library("SpatialExperiment")

####################################### Export out of Seurat
# Are spatial coords in centroids, or in dimension reduction = centroids_exist
Seurat_to_csv <- function(seurat_obj, export_dir){
  
  if (!file.exists(export_dir)){
    dir.create(file.path(export_dir))
  }
  
  # Gene Counts Matrix
  counts_mat <- as.matrix(GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj))) %>% t()
  write.table(counts_mat, paste(export_dir, "/X.csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
  
  #Metadata Matrix
  obs_mat <- seurat_obj[[]]
  write.table(obs_mat, paste(export_dir,"/obs.csv", sep = ""), sep = ",")
  
  # Feature Metadata
  var_mat <- matrix(TRUE, nrow = length(Features(seurat_obj)), ncol = 1)
  rownames(var_mat) <- Features(seurat_obj)
  write.table(var_mat, paste(export_dir,"/var.csv", sep = ""), sep = ",")
  
  # Reductions and Spatial Coordinates
  X_pca <- seurat_obj@reductions$pca@cell.embeddings
  colnames(X_pca) <- sapply(1:ncol(X_pca), function(x){paste("X_pca", x, sep = "")})
  X_umap <- seurat_obj@reductions$umap@cell.embeddings
  colnames(X_umap) <- sapply(1:ncol(X_umap), function(x){paste("X_umap", x, sep = "")})
  
  spatial_list<- lapply(Images(seurat_obj), function(x){seurat_obj@images[[x]]$centroids@coords})
  spatial <- do.call("rbind", spatial_list)
  colnames(spatial) <- sapply(1:ncol(spatial), function(x){paste("spatial",x,sep = "")})
  
  obsm_mat <- cbind(X_pca,X_umap,spatial)
  write.table(obsm_mat,paste(export_dir, "/obsm.csv", sep = ""), row.names = FALSE, sep = ",")}

########################################################
# Import into Seurat
csv_to_Seurat <- function(import_dir, assay_name){
  counts_mat <- read.table(paste(import_dir, "/X.csv", sep = ""), header = FALSE, sep = ",", fill = TRUE) %>% as.matrix() %>% t()
  obs_mat <- read.table(paste(import_dir,"/obs.csv", sep = ""), header = TRUE, sep = ",", row.names = 1, fill = TRUE)
  var_mat <- read.table(paste(import_dir,"/var.csv", sep = ""), header = TRUE, sep = ",", row.names = 1, fill = TRUE)
  obsm_mat <- read.table(paste(import_dir, "/obsm.csv", sep = ""), header = TRUE, sep = ",", fill = TRUE)
  
  rownames(counts_mat) <- rownames(var_mat)
  colnames(counts_mat) <- rownames(obs_mat)
  
  X_pca <- select(obsm_mat, starts_with('X_pca'))
  colnames(X_pca) <- sapply(1:ncol(X_pca), function(x){paste("pca", x, sep = "_")})
  rownames(X_pca) <- rownames(obs_mat)
  
  X_umap <- select(obsm_mat,starts_with('X_umap'))
  colnames(X_umap) <- sapply(1:ncol(X_umap), function(x){paste("umap", x, sep = "_")})
  rownames(X_umap) <- rownames(obs_mat)
  
  spatial <- select(obsm_mat, starts_with("spatial"))
  colnames(spatial) <- sapply(1:ncol(spatial), function(x){paste("spatial", x, sep = "_")})
  rownames(spatial) <- rownames(obs_mat)
  
  seurat_import <- CreateSeuratObject(counts = counts_mat, assay = assay_name)
  seurat_import <- AddMetaData(seurat_import, obs_mat)
  seurat_import[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(X_pca), key = "pca_", assay = DefaultAssay(seurat_import))
  seurat_import[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(X_umap), key = "umap_", assay = DefaultAssay(seurat_import))
  seurat_import[["spatial"]] <- CreateFOV(CreateCentroids(spatial[,c(1,2)]), assay = assay_name)
  
  return(seurat_import)
}

Seurat_to_SPE <- function(seurat_obj){
  spE <- SpatialExperiment(assay = GetAssayData(seurat_obj), colData = seurat_obj[[]], spatialCoords = seurat_obj@images$image$centroids@coords)
  return(spE)
}

Seurat_to_SPIATSPE <- function(seurat_obj){
  spE <- format_image_to_spe(format = "general", 
                             intensity_matrix = as.matrix(GetAssayData(seurat_obj)),
                             phenotypes = as.character(Idents(seurat_obj)), 
                             coord_x = seurat_obj@images$image$centroids@coords[,1],coord_y = seurat_obj@images$image$centroids@coords[,2])
  return(spE)
}
