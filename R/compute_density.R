#' Perform kernel density estimation on SpatialExperiment
#'
#' @param spe A SpatialExperiment object.
#' @param bandwidth The smoothing bandwidth. By default performing automatic bandwidth
#' @param weights Optional weights to be attached to the points.
#' @param kernel The smoothing kernel.Options are gaussian, epanechnikov, quartic or disc.
#' @param scale Scaling factor for density.
#' @param ngrid_x Number of grids in the x-direction.
#' @param ngrid_y Number of grids in the y-direction.
#' @param grid_length_x Grid length in the x-direction.
#' @param grid_length_y Grid length in the y-direction.
#' @param xlim The range of the x-coordinates of the image.
#' @param ylim The range of the y-coordinates of the image.
#' 
#' selection using cross-validation using function spatstat.explore::bw.diggle.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
compute_density <- function(spe, bandwidth = NULL, weights = NULL,
                            kernel = c("gaussian", "epanechnikov", "quartic", "disc"), scale=1e4,
                            ngrid_x = 100, ngrid_y = NULL, 
                            grid_length_x = NULL, grid_length_y = NULL, xlim, ylim
){
  
  if(length(kernel) == 4){
    kernel <- "gaussian"
  }
  
  sc <- SpatialExperiment::spatialCoords(spe)
  
  y <- spatstat.geom::ppp(sc[,1], sc[,2], xlim, ylim)
  
  if(is.null(bandwidth)){
    bandwidth <- spatstat.explore::bw.diggle(y)
  }

  if(!is.null(grid_length_x) | !is.null(grid_length_y)){
    if(is.null(grid_length_y)) 
      grid_length_y <- grid_length_x
    if(is.null(grid_length_x)) 
      grid_length_x <- grid_length_y

    ngrid_x <- round(diff(xlim) / grid_length_x)
    ngrid_y <- round(diff(ylim) / grid_length_y)
  } else {
    if(is.null(ngrid_y))
      ngrid_y <- round(diff(ylim) / diff(xlim) * ngrid_x)
  }

  density_est <- spatstat.explore::density.ppp(y, sigma = bandwidth,
                                               kernel = kernel, weights = weights,
                                               at = "pixels", dimyx = c(ngrid_y, ngrid_x))

  grid_density <- density_est$v * scale
  rownames(grid_density) <- density_est$yrow
  colnames(grid_density) <- density_est$xcol
  grid_density <- setNames(reshape2::melt(grid_density), c('y_grid', 'x_grid', 'density'))

  #grid_density[,c(2,1,3)]
  return(list(grid_density = grid_density[,c(2,1,3)], density_est = density_est))
}
