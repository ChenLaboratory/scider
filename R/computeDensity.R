#' Perform kernel density estimation on SpatialExperiment
#'
#' @param spe A SpatialExperiment object.
#' @param mode Choose either points or pixels. Specifying whether to compute the
#' density at a grid pixel location or at at point.
#' @param kernel The smoothing kernel. Options are gaussian, epanechnikov, quartic or disc.
#' @param bandwidth The smoothing bandwidth. By default performing automatic bandwidth
#' selection using cross-validation using function spatstat.explore::bw.diggle.
#' @param weights Optional weights to be attached to the points.
#' @param scale A munmeric vector to scale the density values. 
#' @param ngrid_x Number of grids in the x-direction. Default to 100.
#' @param ngrid_y Number of grids in the y-direction.
#' @param grid_length_x Grid length in the x-direction.
#' @param grid_length_y Grid length in the y-direction.
#' @param xlim The range of the x-coordinates of the image.
#' @param ylim The range of the y-coordinates of the image.
#'
#' @return Output from spatstat.explore::density.ppp.
#' @export
#'
#' @examples
#' 
#' data("xenium_bc_spe")
#' 
#' dens <- computeDensity(spe)
#' 
computeDensity <- function(spe, mode = "pixels", 
                           kernel = "gaussian", 
                           bandwidth = NULL, weights = NULL, scale=1e4,
                           ngrid_x = 100, ngrid_y = NULL, 
                           grid_length_x = NULL, grid_length_y = NULL, 
                           xlim = NULL, ylim = NULL){

  
  if(!mode %in% c("points","pixels"))
    stop("mode must be either pixels or points.")

  if(!kernel %in% c("gaussian", "epanechnikov", "quartic", "disc"))
    stop("kernel must be one of the followings: gaussian, epanechnikov, quartic or disc.")
  
  sc <- SpatialExperiment::spatialCoords(spe)
  
  if(is.null(xlim))
    xlim <- c(min(sc[,1]), max(sc[,1]))
  
  if(is.null(ylim))
    ylim <- c(min(sc[,2]), max(sc[,2]))

  y <- spatstat.geom::ppp(sc[,1], sc[,2], xlim, ylim)
  
  if(is.null(bandwidth))
    bandwidth <- spatstat.explore::bw.diggle(y)

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
                                               at = mode, dimyx = c(ngrid_y, ngrid_x))
  if(mode == "points"){
    return(density_est)
  } else if (mode == "pixels"){
    
    grid_density <- density_est$v * scale
    rownames(grid_density) <- density_est$yrow
    colnames(grid_density) <- density_est$xcol
    
    #grid_density <- setNames(reshape2::melt(grid_density), c('y_grid', 'x_grid', 'density'))
    
    grid_density <- grid_density |>
      as.data.frame() |>
      rownames_to_column("y_grid") |>
      tidyr::pivot_longer(cols = -y_grid, names_to = "x_grid", values_to = "density")
    
    grid_density <- as.data.frame(sapply(grid_density, as.numeric)) |>
      arrange(x_grid, y_grid)
 
    return(list(grid_density = grid_density[,c(2,1,3)], density_est = density_est))
  }
}
