#' Perform kernel density estimation on SpatialExperiment 
#'
#' @param spe A SpatialExperiment object.
#' @param kernel The smoothing kernel. Options are gaussian, epanechnikov,
#' quartic or disc. ONLY GAUSSIAN IS IMPLEMENTED
#' @param bandwidth The smoothing bandwidth. By default performing automatic
#' bandwidth selection using cross-validation using function
#' spatstat.explore::bw.diggle.
#' @param weights Optional weights to be attached to the points.
#' @param ngrid.x Number of grids in the x-direction. 
#' @param grid.length.x Grid length in the x-direction. 
#' Default to 100 (micron).
#' @param xlim The range of the x-coordinates of the image.
#' @param ylim The range of the y-coordinates of the image.
#' @param diggle Logical. If TRUE, use the Jones-Diggle improved edge
#' @param isVisium Logical. If TRUE, fit hexagonal grids to Visium spots.
#' correction. See spatstat.explore::density.ppp() for details.
#'
#' @return Output from spatstat.explore::density.ppp.
#'
#'
computeDensityHex <- function(spe,
                           kernel = c("gaussian"),
                           bandwidth = NULL, weights = NULL,
                           ngrid.x = NULL,
                           grid.length.x = 100,
                           xlim = NULL, ylim = NULL, diggle = FALSE,
                           isVisium=F) {
  kernel = match.arg(kernel)
  
  sc <- SpatialExperiment::spatialCoords(spe)
  
  if (is.null(xlim)) {
    xlim <- range(sc[,1])
  }
  
  if (is.null(ylim)) {
    ylim <- range(sc[,2])
  }
  
  y <- ppp(sc[, 1], sc[, 2], xlim, ylim)
  
  if (is.null(bandwidth)) {
    bandwidth <- bw.diggle(y)
  }
  
  if (!is.null(grid.length.x)) {
    ngrid.x <- diff(xlim) / grid.length.x
  }
  
  if (isVisium && grid.length.x%%100!=0) {
    warning("For Visium, grid.length.x should be a multiple 
            of 100 to exactly align each spot to a hexagon")
  }
  
  
  density_est = hexDensity::hexDensity(y,
                           bandwidth=bandwidth,
                           weight = weights,
                           xbins=ngrid.x,
                           diggle=diggle,
                           xbnds=xlim,
                           ybnds=ylim
  )
  # Multiply by area 
  density = density_est@count*((diff(density_est@xbnds)/ngrid.x)**2*sqrt(3)/2)
  coords = hexbin::hcell2xy(density_est)
  node_x = (density_est@cell-1)%%density_est@dimen[2]+1
  node_y = (density_est@cell-1)%/%density_est@dimen[2]+1
  
  return(list(
    grid_density = S4Vectors::DataFrame(x_grid=coords$x,
                                        y_grid=coords$y,
                                        node_x=node_x,
                                        node_y=node_y,
                                        density=density),
    density_est=density_est
  ))
}