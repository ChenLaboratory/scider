#' Perform kernel density estimation on SpatialExperiment
#'
#' @param xy A numeric matrix of spatial coordinates.
#' @param kernel The smoothing kernel. Options are gaussian, epanechnikov,
#' quartic or disc.
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
computeDensity <- function(xy,
                           kernel = c("gaussian", "epanechnikov", "quartic", "disc"),
                           bandwidth = NULL, weights = NULL,
                           ngrid.x = NULL,
                           xlim = NULL, ylim = NULL, diggle = FALSE,
                           gridInfo=FALSE) {
  kernel <- match.arg(kernel)
  
  y <- ppp(xy[, 1], xy[, 2], xlim, ylim)
  
  # Using ngrid instead of grid.length.x for backward consistency
  grid.length.x <- diff(xlim)/ngrid.x
  ngrid.x <- round(ngrid.x)
  ngrid.y <- round(diff(ylim)/grid.length.x)
  
  density_est <- density.ppp(y,
                             sigma = bandwidth,
                             kernel = kernel, weights = weights,
                             dimyx = c(ngrid.y, ngrid.x),
                             diggle = diggle
  )
  grid_density <- density_est$v * density_est$xstep * density_est$ystep
  # Reshape data from wide to long
  rownames(grid_density) <- density_est$yrow
  colnames(grid_density) <- density_est$xcol
  grid_density <- as.data.frame(grid_density) |>
    rownames2col("y_grid")
  
  reshaped_dat <- stats::reshape(grid_density,
                                 direction = "long",
                                 varying = names(grid_density)[(names(grid_density) != "y_grid")],
                                 v.names = "density",
                                 timevar = "x_grid",
                                 times = names(grid_density)[(names(grid_density) != "y_grid")]
  )
  
  reshaped_dat <- reshaped_dat[, names(reshaped_dat) != "id"]
  rownames(reshaped_dat) <- NULL
  grid_density <- as.data.frame(vapply(
    reshaped_dat, as.numeric,
    numeric(nrow(reshaped_dat))
  ))
  grid_density <- grid_density[order(
    grid_density$x_grid,
    grid_density$y_grid
  ), ]
  
  if (gridInfo) {
    ngrid.x <- density_est$dim[2]
    ngrid.y <- density_est$dim[1]
    grid_density = grid_density[, c(2, 1)]
    # horizontal ind
    grid_density$node_x <- rep(seq_len(ngrid.x),each = ngrid.y)
    # vertical ind
    grid_density$node_y <- rep(seq_len(ngrid.y),ngrid.x)
    grid_density$node <- paste(grid_density$node_x,
                               grid_density$node_y,
                               sep = "-")
    return(list(
      grid_density = S4Vectors::DataFrame(grid_density),
      grid_info = list(
        dims = c(ngrid.x, ngrid.y),
        xlim = xlim,
        ylim = ylim,
        xcol = density_est$xcol,
        yrow = density_est$yrow,
        xstep = density_est$xstep,
        ystep = density_est$ystep,
        bandwidth = bandwidth,
        grid_type = "square"
      )
    ))
  } else {
    # Multiply by area
    return(grid_density$density)
  }
}
