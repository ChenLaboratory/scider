#' Plot grid-based density.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' @param probs Numeric value between 0 and 1, used for filtering uninformative grid, default is 0.8.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#' 
#' spe <- gridDensity(spe)
#'
#' plotDensity(spe, coi = "Breast cancer")
#' 
#' plotDensity(spe, coi = "Fibroblasts")
#' 

plotDensity <- function(spe, coi, probs = 0.8){

  grid_data <- spe@metadata$grid_density

  coi_clean <- janitor::make_clean_names(coi)
  dens_cols <- paste("density", coi_clean, sep="_")
    
  if(!all(dens_cols %in% colnames(grid_data)))
    stop("Density of COI is not yet computed.")

  grid_data$density_coi_average <- rowMeans(grid_data[, which(colnames(grid_data) %in% dens_cols), drop=FALSE])

  kp <- grid_data$density_coi_average >= quantile(grid_data$density_coi_average, probs=probs)

  p <- ggplot(grid_data[kp,], aes(x=x_grid, y=y_grid, z = density_coi_average)) +
    geom_raster(aes(fill = density_coi_average)) + 
    theme_classic() +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral"))) + 
    labs(x = "x", y = "y", fill = "Density") +
    ggtitle(coi) 

  return(p)
}