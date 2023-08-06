#' Perform kernel density estimation on SpatialExperiment for cell types of interest
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs). 
#' @param bandwidth Smoothing parameter
#' @param coi.cn A character vector. Column name of the cell types. 
#' By default is cell_type.
#' @param ngrid_x Number of grids in the x-direction.
#' @param ngrid_y Number of grids in the y-direction.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#' 
#' spe_dens <- grid_dens_by_coi(spe, coi = c("Breast cancer", "Fibroblasts"))
#' 
#'
grid_dens_by_coi <- function(spe, 
                          coi, 
                          coi.cn = "cell_type",
                          ngrid_x = 100, ngrid_y = 100, 
                          bandwidth = 50) {
  
  # define canvas
  spatialCoordsNames(spe) <- c("x_centroid", "y_centroid")
  coord <- spatialCoords(spe)
  xlim <- c(min(coord[, "x_centroid"]), max(coord[, "x_centroid"]))
  ylim <- c(min(coord[, "y_centroid"]), max(coord[, "y_centroid"]))
  
  coldat <- as.data.frame(colData(spe))
  
  # compute density for each cell type and then, filter
  dens_list <- list()
  for(ii in 1:length(coi)){
    
    # subset data to this COI
    obj <- spe[, coldat[,coi.cn] == coi[ii]]
    
    # compute density
    spe_den <- compute_density(obj, 
                           bandwidth = bandwidth, 
                           mode = "pixels", 
                           xlim = xlim, 
                           ylim = ylim, 
                           ngrid_x = ngrid_x,
                           ngrid_y = ngrid_y)
    
    den_i <- metadata(spe_den)$grid_density
    
    colnames(den_i)[colnames(den_i) == "density"] <- paste0(gsub(" ", "_",coi[ii]),"_density")
    
    # give grid names
    den_i <- den_i %>%
      mutate(node_x = rep(1:ngrid_x, ngrid_y),# row ind
             node_y = rep(1:ngrid_y, each = ngrid_x),# col ind
             node = paste0(node_x,"-",node_y)) %>%
      dplyr::select(c("x_grid","y_grid","node_x","node_y","node",everything()))
    
    dens_list[[ii]] <- den_i
  }
  
  dens_result <- dens_list %>%
    purrr::reduce(left_join, by = c("x_grid","y_grid","node_x","node_y","node"))
  
  # return
  
  metadata(spe)$grid_density <- dens_result
  metadata(spe)$grid_info <- metadata(spe_den)$grid_info
  
  return(spe)
  
}