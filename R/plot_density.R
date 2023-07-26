#' Plot grid-based density.
#'
#' @param spe A SpatialExperiment object.
#' @param q.tile Numeric value between 0 and 1, used for filtering uninformative grid, default is 0.8.
#'  
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#' 
#' spe <- compute_density(spe)
#'
#' plot_density(spe)
plot_density <- function(spe, q.tile = 0.8){
  
  den <- S4Vectors::metadata(spe)$grid_density
  
  filter <- den$density > quantile(den$density, q.tile)
  
  p <- ggplot(den[filter,], aes(x=x_grid, y=y_grid, z = density)) +
    geom_raster(aes(fill = density)) + 
    theme_classic() +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral")))
  
  return(p)
}






