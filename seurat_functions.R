library("Seurat")
library("readr")
library("dplyr")
library("stringr")

####################################### Export out of Seurat
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
  if ("pca" %in% names(seurat_obj)){
  X_pca <- seurat_obj@reductions$pca@cell.embeddings
  colnames(X_pca) <- sapply(1:ncol(X_pca), function(x){paste("X_pca", x, sep = "")})}
  if ("umap" %in% names(seurat_obj)){
  X_umap <- seurat_obj@reductions$umap@cell.embeddings
  colnames(X_umap) <- sapply(1:ncol(X_umap), function(x){paste("X_umap", x, sep = "")})}
  
  if (!is.null(Images(seurat_obj))){
  spatial_list<- lapply(Images(seurat_obj), function(x){seurat_obj@images[[x]]$centroids@coords})
  spatial <- do.call("rbind", spatial_list)
  spatial <- spatial[,c(2,1)]
  colnames(spatial) <- sapply(1:ncol(spatial), function(x){paste("spatial",x,sep = "")})
  
  if (.hasSlot(seurat_obj@images[[DefaultFOV(seurat_obj)]], "image")){
    if (!file.exists(paste(export_dir,"/uns",sep = ''))){
      dir.create(file.path(paste(export_dir,'/uns', sep = '')))
    }
    
    spatial_r <- as.matrix(seurat_obj@images[[DefaultFOV(seurat_obj)]]@image[,,1]) * 255
    spatial_g <- as.matrix(seurat_obj@images[[DefaultFOV(seurat_obj)]]@image[,,2]) * 255
    spatial_b <- as.matrix(seurat_obj@images[[DefaultFOV(seurat_obj)]]@image[,,3]) * 255
    write.table(spatial_r,paste(export_dir, "/uns/spatial_lowres_r.csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
    write.table(spatial_g,paste(export_dir, "/uns/spatial_lowres_g.csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
    write.table(spatial_b,paste(export_dir, "/uns/spatial_lowres_b.csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
    
    scale_factor_list <- as.list(ScaleFactors(seurat_obj@images[[DefaultFOV(seurat_obj)]]))
    class(scale_factor_list) <- NULL
    name_map <- c(spot = 'spot_diameter_fullres', fiducial = 'fiducial_diameter_fullres', hires = 'tissue_hires_scalef', lowres = 'tissue_lowres_scalef')
    names(scale_factor_list) <- recode(names(scale_factor_list), !!!name_map)
    write.table(data.frame(scale_factor_list), paste(export_dir, "/uns/spot_size.csv", sep = ""), row.names = FALSE, sep = ",")
  }
  }
  
  vars <- c("X_pca", "X_umap", "spatial")
  # Get only the ones that exist
  env <- environment()
  existing_vars <- vars[sapply(vars, exists, envir = env)]
  # Combine them if any exist
  if (length(existing_vars) > 0) {
    obsm_mat <- do.call(cbind, mget(existing_vars))
    write.table(obsm_mat,paste(export_dir, "/obsm.csv", sep = ""), row.names = FALSE, sep = ",")
  }
}

########################################################
# Import into Seurat
csv_to_Seurat <- function(import_dir, assay_name){
  counts_mat <- read.table(paste(import_dir, "/X.csv", sep = ""), header = FALSE, sep = ",", fill = TRUE) %>% as.matrix() %>% t()
  obs_mat <- read.table(paste(import_dir,"/obs.csv", sep = ""), header = TRUE, sep = ",", row.names = 1, fill = TRUE)
  var_mat <- read.table(paste(import_dir,"/var.csv", sep = ""), header = TRUE, sep = ",", row.names = NULL, fill = TRUE)
  rownames(var_mat) <- make.unique(var_mat[,1])
  var_mat <- var_mat[,-1]
  if (file.exists(paste(import_dir,"/obsm.csv", sep = ""))){
    content <- readLines(paste(import_dir, '/obsm.csv', sep = ""), n = 1, warn = FALSE)
    if (content != ""){
      obsm_mat <- read.table(paste(import_dir, "/obsm.csv", sep = ""), header = TRUE, sep = ",", fill = TRUE)}
  }
  
  if (file.exists(paste(import_dir, "/uns/spot_size.csv", sep = ''))){
    
    scale_factors <- read.table(paste(import_dir,'/uns/spot_size.csv', sep = ''), sep = ',', header = TRUE)
    
    if (file.exists(paste(import_dir, "/uns/spatial_lowres_r.csv", sep = ''))){
      spatial_r <- as.array(as.matrix(read.table(paste(import_dir,'/uns/spatial_lowres_r.csv', sep = ''), sep = ',', header = FALSE))) / 255
      spatial_g <- as.array(as.matrix(read.table(paste(import_dir,'/uns/spatial_lowres_g.csv', sep = ''), sep = ',', header = FALSE))) / 255
      spatial_b <- as.array(as.matrix(read.table(paste(import_dir,'/uns/spatial_lowres_b.csv', sep = ''), sep = ',', header = FALSE))) / 255
      
      dimnames(spatial_r) <- NULL
      dimnames(spatial_g) <- NULL
      dimnames(spatial_b) <- NULL
    }
    else if (file.exists(paste(import_dir, "/uns/spatial_hires_r.csv", sep = ''))){
      warning("Low resolution image missing. Using high resolution image instead.")
      
      spatial_r <- as.array(as.matrix(read.table(paste(import_dir,'/uns/spatial_hires_r.csv', sep = ''), sep = ',', header = FALSE))) / 255
      spatial_g <- as.array(as.matrix(read.table(paste(import_dir,'/uns/spatial_hires_g.csv', sep = ''), sep = ',', header = FALSE))) / 255
      spatial_b <- as.array(as.matrix(read.table(paste(import_dir,'/uns/spatial_hires_b.csv', sep = ''), sep = ',', header = FALSE))) / 255
      
      dimnames(spatial_r) <- NULL
      dimnames(spatial_g) <- NULL
      dimnames(spatial_b) <- NULL
    }
    spatial_image <- array(c(spatial_r, spatial_g, spatial_b), dim = c(nrow(spatial_r), ncol(spatial_r), 3))
  }
  
  rownames(counts_mat) <- rownames(var_mat)
  colnames(counts_mat) <- rownames(obs_mat)
  
  seurat_import <- CreateSeuratObject(counts = counts_mat, assay = assay_name)
  seurat_import <- AddMetaData(seurat_import, obs_mat)
  
  if (exists("obsm_mat")){
  if (any(startsWith(names(obsm_mat), "X_pca"))){
  X_pca <- select(obsm_mat, starts_with('X_pca'))
  colnames(X_pca) <- sapply(1:ncol(X_pca), function(x){paste("pca", x, sep = "_")})
  rownames(X_pca) <- rownames(obs_mat)
  seurat_import[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(X_pca), key = "pca_", assay = DefaultAssay(seurat_import))
  }
  
  if (any(startsWith(names(obsm_mat), "X_umap"))){
  X_umap <- select(obsm_mat,starts_with('X_umap'))
  colnames(X_umap) <- sapply(1:ncol(X_umap), function(x){paste("umap", x, sep = "_")})
  rownames(X_umap) <- rownames(obs_mat)
  seurat_import[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(X_umap), key = "umap_", assay = DefaultAssay(seurat_import))
  }
  
  if (any(startsWith(names(obsm_mat), "spatial"))){
  spatial <- select(obsm_mat, starts_with("spatial"))
  spatial[,c(1,2)] <- spatial[,c(2,1)]
  colnames(spatial) <- sapply(1:ncol(spatial), function(x){paste("spatial", x, sep = "_")})
  rownames(spatial) <- rownames(obs_mat)
  fov_curr <- CreateFOV(CreateCentroids(spatial[,c(1,2)]), assay = assay_name)
  
  if (exists('spatial_image')){
    visium_new = new("VisiumV2", image = spatial_image, scale.factors = scalefactors(spot = scale_factors$spot_diameter_fullres, fiducial = scale_factors$fiducial_diameter_fullres, hires = scale_factors$tissue_hires_scalef, lowres = scale_factors$tissue_lowres_scalef), boundaries = list(centroids = CreateCentroids(spatial[,c(1,2)])), assay = assay_name, key = Key(fov_curr))
    seurat_import[['spatial']] <- visium_new
  }
  else{
    seurat_import[["spatial"]] <- fov_curr
  }
  }
  }
  
  return(seurat_import)
}

Seurat_to_SPE <- function(seurat_obj){
  library("SpatialExperiment")
  spE <- SpatialExperiment(assay = GetAssayData(seurat_obj), colData = seurat_obj[[]], spatialCoords = seurat_obj@images$image$centroids@coords)
  return(spE)
}

Seurat_to_SPIATSPE <- function(seurat_obj){
  library('SPIAT')
  spE <- format_image_to_spe(format = "general", 
                             intensity_matrix = as.matrix(GetAssayData(seurat_obj)),
                             phenotypes = as.character(Idents(seurat_obj)), 
                             coord_x = seurat_obj@images$image$centroids@coords[,1],coord_y = seurat_obj@images$image$centroids@coords[,2])
  return(spE)
}
