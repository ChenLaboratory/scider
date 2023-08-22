#' Plot contour lines.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' @param ... Aesthetic mappings to pass to `ggplot2::aes_string()`.
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
#' coi <- "Breast cancer"
#' 
#' spe <- getContour(spe, coi = coi)
#' 
#' plotContour(spe, coi = coi, size = 0.3, alpha = 0.2)
#' 

plotContour <- function(spe, coi, ...){

  coi_clean <- janitor::make_clean_names(coi)
  coi_clean_contour <- paste(coi_clean, "contour", sep = "_")

  contour_data <- spe@metadata[[coi_clean_contour]]
  levs <- unique(contour_data$level)
  nlevs <- length(levs)

  p <- plotSpatial(spe, ...) + 
    ggplot2::geom_path(data = contour_data, 
      ggplot2::aes(x = x, y = y, group = group, color = level_factor)) + 
    scale_color_manual(values= rev(RColorBrewer::brewer.pal(nlevs, "Spectral"))) +
    theme_classic() +
    labs(x = "x", y = "y", color = "Density") +
    ggtitle(coi) 

  return(p)
}

