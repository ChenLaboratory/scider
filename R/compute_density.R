#' Perform kernel density estimation on SpatialExperiment
#'
#' @param spe A SpatialExperiment object.
#' @param bandwidth The smoothing bandwidth. By default performing automatic bandwidth
#' @param weights Optional weights to be attached to the points.
#' @param kernel The smoothing kernel.Options are gaussian, epanechnikov, quartic or disc.
#' @param den_scaler Numeric value, used for scaling density value, default is 1e4.
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
#' 
#' data("xenium_bc_spe")
#' 
#' compute_density(spe)
#' 
compute_density <- function(spe, mode = c("points","pixels"), bandwidth = NULL, weights = NULL,
                            kernel = c("gaussian", "epanechnikov", "quartic", "disc"), 
                            den_scaler = 1e4, ngrid_x = 100, ngrid_y = 100, 
                            grid_length_x = NULL, grid_length_y = NULL,
                            xlim = NULL, ylim = NULL){

  if(length(kernel) == 4){
    kernel <- "gaussian"
  }
  
  if(length(mode) == 2){
    mode <- "pixels"
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
    
    colData(spe)$density = density_est
    
  } else if (mode == "pixels"){
    
    grid_density <- density_est$v
    
    rownames(grid_density) <- density_est$yrow
    
    colnames(grid_density) <- density_est$xcol
    
    grid_density <- grid_density |> 
      as.data.frame() |> 
      tibble::rownames_to_column("y_grid") |> 
      tidyr::pivot_longer(-y_grid, names_to = "x_grid", values_to = "density") |> 
      dplyr::select(c("x_grid","y_grid","density")) |>
      dplyr::mutate(density = density * den_scaler,
                    x_grid = as.numeric(x_grid),
                    y_grid = as.numeric(y_grid)) |>
      as.data.frame()
    
    rownames(grid_density) <- paste0("grid_", seq(nrow(grid_density)))

    metadata(spe) <- list("grid_density" = grid_density)
    
  }

  return(spe)
}
