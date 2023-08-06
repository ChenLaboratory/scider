#' Plot ROIs on spatial.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs). 
#' @param col.cn A character vector. Column name of the cell types. 
#' By default is cell_type.
#' @param n Integer value. The number of distinct color to be generated, default is 30.
#' @param threshold Integer value. The threshold (number of grid) used to find large ROIs.
#' @param showlegend Logical. Show legend or not.
#' @param ... Aesthetic mappings to pass to `ggplot2::aes_string()` for point.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' 
#' data("xenium_bc_spe")
#' 
#' spe_dens <- grid_dens_by_coi(spe, coi = c("Breast cancer", "Fibroblasts"))
#' 
#' spe_rois <- find_roi_by_dens(spe_dens, coi = c("Breast cancer", "Fibroblasts"), 
#' method = "walktrap", steps = 5)
#' 
#' plot_roi(spe_rois, coi = c("Breast cancer", "Fibroblasts"), size = 0.01, alpha = 0.2)
#' 
plot_roi <- function(spe, coi = NULL, col.cn = "cell_type", n = 30, 
                     threshold = 0, showlegend = FALSE, ...){
  rois <- metadata(spe)$components
  
  dat <- as.data.frame(colData(spe))
  
  if (!is.null(coi)) {
    spe <- spe[, dat[,col.cn] %in% coi]
  }
  
  posdat <- as.data.frame(spatialCoords(spe))
  
  dat <- as.data.frame(colData(spe)) %>%
    cbind(posdat)
  
  set.seed(100)
  col.p <- randomcoloR::distinctColorPalette(n)
  
  xlim <- c(min(posdat[,1]), max(posdat[,1]))
  ylim <- c(min(posdat[,2]), max(posdat[,2]))
  plot.xlim <- xlim + c(-1e-10, 1e-10)
  plot.ylim <- ylim + c(-1e-10, 1e-10)
  
  filtered <- which(table(rois$component) > threshold)
  rois_filtered <- rois[rois$component %in% filtered, ]
  
  roi_plot <- plot_spatial(spe, ...) +
    geom_tile(data = rois_filtered, aes(x = xcoord, y = ycoord, fill = component), alpha = 0.6) +
    scale_fill_manual(values = rep(col.p[1:10], 100)) +
    scale_x_continuous(limits = plot.xlim) +
    scale_y_continuous(limits = plot.ylim)
  
  if (isFALSE(showlegend)){
    roi_plot <- roi_plot +
      theme(legend.position = "none")
  } 
  
  return(roi_plot)
}




