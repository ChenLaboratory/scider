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
#' @param ngrid.x Number of grids in the x-direction. Default to 100.
#' @param ngrid.y Number of grids in the y-direction.
#' @param grid.length.x Grid length in the x-direction.
#' @param grid.length.y Grid length in the y-direction.
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
                           ngrid.x = 100, ngrid.y = NULL, 
                           grid.length.x = NULL, grid.length.y = NULL, 
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

  if(!is.null(grid.length.x) | !is.null(grid.length.y)){
    if(is.null(grid.length.y)) 
      grid.length.y <- grid.length.x
    if(is.null(grid.length.x)) 
      grid.length.x <- grid.length.y
    ngrid.x <- round(diff(xlim) / grid.length.x)
    ngrid.y <- round(diff(ylim) / grid.length.y)
  } else {
    if(is.null(ngrid.y))
      ngrid.y <- round(diff(ylim) / diff(xlim) * ngrid.x)
  }

  density_est <- spatstat.explore::density.ppp(y, sigma = bandwidth,
                                               kernel = kernel, weights = weights,
                                               at = mode, dimyx = c(ngrid.y, ngrid.x))
  if(mode == "points"){
    return(density_est)
  } else if (mode == "pixels"){
    
    grid_density <- density_est$v * scale
    rownames(grid_density) <- density_est$yrow
    colnames(grid_density) <- density_est$xcol
    
    grid_density <- as.data.frame(grid_density) |>
      rownames2col("y_grid")
    
    reshaped_dat <- stats::reshape(grid_density, direction = "long", 
                                   varying = names(grid_density)[(names(grid_density) != "y_grid")],
                                   v.names = "density", 
                                   timevar = "x_grid", 
                                   times = names(grid_density)[(names(grid_density) != "y_grid")])
    
    reshaped_dat <- reshaped_dat[,names(reshaped_dat) != "id"]
    rownames(reshaped_dat) <- NULL
    grid_density <- as.data.frame(sapply(reshaped_dat, as.numeric))
    grid_density <- grid_density[order(grid_density$x_grid, grid_density$y_grid),]
      
 
    return(list(grid_density = S4Vectors::DataFrame(grid_density[,c(2,1,3)]), density_est = density_est))
  }
}
