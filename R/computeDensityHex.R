#' Perform kernel density estimation on SpatialExperiment 
#'
#' @param xy A numeric matrix of spatial coordinates.
#' @param kernel The smoothing kernel. Options are gaussian, epanechnikov,
#' quartic or disc. ONLY GAUSSIAN IS IMPLEMENTED
#' @param bandwidth The smoothing bandwidth. By default performing automatic
#' bandwidth selection using cross-validation using function
#' spatstat.explore::bw.diggle.
#' @param weights Optional weights to be attached to the points.
#' @param ngrid.x Number of grids in the x-direction. 
#' @param xlim The range of the x-coordinates of the image.
#' @param ylim The range of the y-coordinates of the image.
#' @param diggle Logical. If TRUE, use the Jones-Diggle improved edge
#' correction. See spatstat.explore::density.ppp() for details.
#' @param gridInfo Logical. If TRUE, then the grid information is also returned.
#'
#' @return Output from spatstat.explore::density.ppp.
#'
#'
computeDensityHex <- function(xy,
                              kernel = c("gaussian"),
                              bandwidth = NULL, weights = NULL,
                              ngrid.x = NULL,
                              xlim = NULL, ylim = NULL, diggle = FALSE,
                              gridInfo=FALSE) {
  kernel <- match.arg(kernel)
  
  # explicitly providing x,y to hexDensity is faster.
  density_est <- hexDensity::hexDensity(x=as.double(xy[,1]),
                                        y=as.double(xy[,2]),
                                        bandwidth = bandwidth,
                                        weight = weights,
                                        xbins = ngrid.x,
                                        edge = !gridInfo,
                                        diggle = diggle,
                                        xbnds = xlim,
                                        ybnds = ylim)
  if (gridInfo) {
    coords <- hexbin::hcell2xy(density_est)
    node_x <- (density_est@cell-1)%%density_est@dimen[2]+1
    node_y <- (density_est@cell-1)%/%density_est@dimen[2]+1
    # Max of node_y might be 1 or 2 less than dimen[1] which is ok for 
    # hexbin/hexDensity but can be a bit unintuitive to work with. 
    density_est@dimen[1] <- max(node_y) 
    return(list(
      grid_density = S4Vectors::DataFrame(x_grid=coords$x,
                                          y_grid=coords$y,
                                          node_x=node_x,
                                          node_y=node_y,
                                          node=paste(node_x,node_y,sep="-")),
      grid_info = list(
        dims = density_est@dimen[2:1],
        xlim = xlim,
        ylim = ylim,
        xstep = diff(xlim)/density_est@xbins,
        ystep = (diff(ylim)*sqrt(3))/(2*density_est@shape*density_est@xbins),
        xbins = density_est@xbins,
        shape = density_est@shape,
        bandwidth = bandwidth,
        grid_type = "hex"
      )
    ))
  } else {
    # Multiply by area
    return(density_est@count*((diff(density_est@xbnds)/ngrid.x)**2*sqrt(3)/2))
  }
}
