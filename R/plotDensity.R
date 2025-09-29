#' Plot grid-based density.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs) to be plotted.
#' Default to all cell types.
#' @param probs Numeric value between 0 and 1, used for filtering
#' uninformative grid, default is 0.5.
#' @param ... Parameters pass to \link[scider]{plotGrid}
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
plotDensity <- function(spe, coi = NULL, probs = 0.5,...) {
  coi_clean <- `if`(is.null(coi),"overall",cleanName(coi))
  dens_cols <- paste("density", coi_clean, sep = "_")
  
  if (!all(dens_cols %in% colnames(spe@metadata$grid_density))) {
    stop("Density of COI is not yet computed.")
  }
  spe@metadata$grid_density$density_coi_average <- rowSums(as.matrix(
    spe@metadata$grid_density[, dens_cols, drop = FALSE]
  ))
  plotGrid(spe,
           group.by="density_coi_average",
           probs=probs,
           label="Density",
           ...) +
    ggtitle(paste(coi_clean, collapse=", "))
}

utils::globalVariables(c("x_grid", "y_grid", "density_coi_average"))