#' Perform kernel density estimation on SpatialExperiment for cell types of interest
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs). Default to all cell types.
#' @param id A character. The name of the column of colData(spe) containing the cell type identifiers.
#' Set to cell_type by default.
#' @param kernel The smoothing kernel. Options are gaussian, epanechnikov, quartic or disc.
#' @param bandwidth The smoothing bandwidth. By default performing automatic bandwidth
#' selection using cross-validation using function spatstat.explore::bw.diggle.
#' @param scale A munmeric vector to scale the density values. 
#' @param ngrid_x Number of grids in the x-direction. Default to 100.
#' @param ngrid_y Number of grids in the y-direction. 
#' @param grid_length_x Grid length in the x-direction.
#' @param grid_length_y Grid length in the y-direction.
#'
#' @return A SpatialExperiment object. Grid density estimates for all cell type of interest are stored in metadata(spe)$grid_density. Grid information is stored in metadata(spe)$grid_info
#' 
#' @export
#'
#' @examples
#' 
#' data("xenium_bc_spe")
#' 
#' gridDensity(spe)
#'

gridDensity <- function(spe, 
                        coi = "all", 
                        id = "cell_type",
                        kernel = "gaussian",
                        bandwidth = NULL,
                        scale = 1e4, 
                        ngrid_x = 100, ngrid_y = NULL, 
                        grid_length_x = NULL, grid_length_y = NULL
                        ) {
  
  if(! id %in% colnames(colData(spe))) 
    stop(paste(id, "is not a column of the colData."))

  if(coi == "all") 
    coi <- names(table(colData(spe)[[id]]))

  if(any(! coi %in% names(table(colData(spe)[[id]]))))
    stop("One or more cell types of interest is not found.")

  coi_clean <- janitor::make_clean_names(coi)
  
  # define canvas
  spatialCoordsNames(spe) <- c("x_centroid", "y_centroid")
  coord <- spatialCoords(spe)
  xlim <- c(min(coord[, "x_centroid"]), max(coord[, "x_centroid"]))
  ylim <- c(min(coord[, "y_centroid"]), max(coord[, "y_centroid"]))
  
  if(is.null(metadata(spe))) metadata(spe) <- list()

  # Reset when the function is rerun again
  metadata(spe)$grid_density <- metadata(spe)$grid_info <- NULL

  # compute density for each cell type and then, filter
  for(ii in 1:length(coi)){
    
    # subset data to this COI
    sub <- grep(coi[ii], colData(spe)[[id]])
    obj <- spe[, sub]
    
    # compute density
    out <- computeDensity(obj, mode = "pixels", kernel = kernel, 
                           bandwidth = bandwidth, scale = scale, 
                           ngrid_x = ngrid_x, ngrid_y = ngrid_y, 
                           grid_length_x = grid_length_x, grid_length_y = grid_length_y,
                           xlim = xlim, ylim = ylim)
    RES <- out$grid_density
    
    ngrid_x <- out$density_est$dim[2]
    ngrid_y <- out$density_est$dim[1]
    
    if(is.null(metadata(spe)$grid_density)){
      metadata(spe) <- list("grid_density" = RES[, 1:2])
      metadata(spe)$grid_density <- RES[, 1:2]
      metadata(spe)$grid_density$node_x <- rep(1:ngrid_x, each = ngrid_y) # horizontal ind
      metadata(spe)$grid_density$node_y <- rep(1:ngrid_y, ngrid_x) # vertical ind
      metadata(spe)$grid_density$node <- paste(metadata(spe)$grid_density$node_x,
                                               metadata(spe)$grid_density$node_y, sep = "-")
    }
    metadata(spe)$grid_density <- cbind(metadata(spe)$grid_density, RES$density)
    colnames(metadata(spe)$grid_density)[5 + ii] <- paste("density", coi_clean[ii], sep="_")
    
    # grid info
    if(is.null(metadata(spe)$grid_info)){
      spe@metadata$grid_info <- list(dims = c(ngrid_x, ngrid_y), 
                                     xlim = xlim, 
                                     ylim = ylim, 
                                     xcol = out$density_est$xcol, 
                                     yrow = out$density_est$yrow, 
                                     xstep = out$density_est$xstep, 
                                     ystep = out$density_est$ystep)
    }
  }
  return(spe)
}
