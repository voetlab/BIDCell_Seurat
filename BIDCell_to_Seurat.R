library(ggplot2)
library(dplyr)
library(sp)
library(data.table)
library(Seurat)
library(patchwork)
library(raster)
library(terra)
library(magrittr)
library(sf)
options(Seurat.object.assay.version = "v5")

#run BIDcell beforehand https://github.com/SydneyBioX/BIDCell
#MWork in Progress
#Made by Markus Boesch and Francisco Pestana
#name=output name of file and FOV
#BIDCell_outs_folder=location of the BIDCell output
#output_folder=location where Seurat rds is suppose to be saved

BIDCell_to_Seurat<-function(name,BIDCell_outs_folder, output_folder){
  folder<-list.dirs(paste0(BIDCell_outs_folder,"/model_outputs", sep=""))[2] #get the last model folder       
    rasterLayer <- flip(raster(paste0(folder,"/test_output/epoch_1_step_4000_connected.tif", sep="")),"y") #import tif to get the segmentation
  polys <- rasterToPolygons(rasterLayer, fun=function(x){x>0}, dissolve=TRUE,n = 4) #generate Polygons 
  sf_polys <- st_as_sf(polys,point = FALSE, merge = T, connect8 = T)
  sf_polys <- sf_polys[order(sf_polys$epoch_1_step_4000_connected),]
  sf_polys$cellID <- as.character(sf_polys$epoch_1_step_4000_connected)
  # Convert sf object to a dataframe and include geometry details
  df <- sf_polys %>%
    st_coordinates() %>%
    as.data.frame()
  # # Create dataframe
  df_boundaries <- data.frame(
    cell_id = df$L3,
    vertex_x = df$X, # Assuming the first column represents the X coordinate
    vertex_y = df$Y  # Assuming the second column represents the Y coordinate
  )
  df_boundaries <- df_boundaries[order(df_boundaries$cell_id),]
  # Assuming your dataframe is named df
  df_boundaries <- df_boundaries %>%
    group_by(cell_id) %>%
    distinct(vertex_x, vertex_y, .keep_all = TRUE) %>%
    ungroup()
  
  #get Centroids info
  centroids<-st_centroid(sf_polys)
  df2<-centroids %>%
    st_coordinates() %>%
    as.data.frame()
  df_centroids<-data.frame(cell=rownames(df2), #because IDs are Numbers
                           x= df2$X, # Assuming the first column represents the X coordinate
                           y = df2$Y 
  )
  expr_mat<-read.csv(paste0(BIDCell_outs_folder,"/cell_gene_matrices/",stringr::str_sub(folder, start= -19),"/expr_mat.csv", sep=""))
  molecules<-read.csv(paste0(BIDCell_outs_folder,"/transcripts.csv.gz", sep=""))
  
  segmentations.data <- list( "centroids"=CreateCentroids(df_centroids) ,"segmentation" = CreateSegmentation(df_boundaries))
  df_molecules<-data.frame("x"=molecules$x_location, "y"=molecules$y_location, "gene"=molecules$feature_name) 
  coords <- CreateFOV( #generate a FOV based on BIDCell polygons
    coords = segmentations.data,
    type = c("centroids", "segmentation"),
    molecules = df_molecules,
    assay = "Xenium"
  )
  expr_mat$X<-NULL
  expr_mat<-t(expr_mat)
  colnames(expr_mat)<-expr_mat[1,]
  expr_mat<-expr_mat[-1,]
  expr_mat <- as(expr_mat, "dgCMatrix")
  xenium.obj = CreateSeuratObject(counts = expr_mat, assay="Xenium") #generate SeuratObject
  xenium.obj[[name]] <- coords
  xenium.obj<-subset(xenium.obj, subset=nCount_Xenium >10) #minimum QC
  xenium.obj$technique<-"Xenium"
  saveRDS(xenium.obj, file = paste0(output_folder,"/",name,'.rds', sep=""))
}
