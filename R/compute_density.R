#' Perform kernel density estimation on SpatialExperiment
#'
#' @param spe A SpatialExperiment object.
#' @param mode Choose either points or pixels. Specifying whether to compute the
#' density at a grid pixel location or at at point.
#' @param grid_size Grid size to use when the mode is set to pixels.
#' @param bandwidth The smoothing bandwidth. By default performing automatic bandwidth
#' selection using cross-validation using function spatstat.explore::bw.diggle.
#' @param weights Optional weights to be attached to the points.
#' @param kernel The smoothing kernel.Options are gaussian, epanechnikov, quartic or disc.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#' 
#' data("xenium_bc_spe")
#' 
#' compute_density(spe)
#' 
compute_density <- function(spe, mode = c("points","pixels"), bandwidth = NULL, weights = NULL,
                            kernel = c("gaussian", "epanechnikov", "quartic", "disc"),
                            grid_size = 100, xlim = NULL, ylim = NULL){

  if(length(kernel) == 4){
    kernel <- "gaussian"
  }

  sc <- SpatialExperiment::spatialCoords(spe)
  
  if(is.null(xlim)){
    xlim <- c(min(sc[,1]), max(sc[,1]))
  }
  
  if(is.null(ylim)){
    ylim <- c(min(sc[,2]), max(sc[,2]))
  }

  y <- spatstat.geom::ppp(sc[,1], sc[,2], xlim, ylim)

  if(is.null(bandwidth)){
    bandwidth <- spatstat.explore::bw.diggle(y)
  }
  
  if(length(mode) == 2){
    mode <- "pixels"
  }

  density_est <- spatstat.explore::density.ppp(y, sigma = bandwidth,
                                               kernel = kernel, weights = weights,
                                               at = mode, dimyx = c(grid_size, grid_size))

  if(mode == "points"){
    colData(spe)$density = density_est

  } else if (mode == "pixels"){
    grid_density <- density_est$v
    rownames(grid_density) <- density_est$yrow
    colnames(grid_density) <- density_est$xcol
    grid_density <- grid_density |> 
      as.data.frame() |> 
      rownames_to_column("y_grid") |> 
      pivot_longer(-y_grid, names_to = "x_grid", values_to = "density") |> 
      dplyr::select(c("x_grid","y_grid","density"))

    metadata(spe) <- list("grid_density" = grid_density)
  }

  return(spe)
}

