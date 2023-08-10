#' Plot ROIs on spatial.
#'
#' @param spe A SpatialExperiment object.
#' @param id A character. The name of the column of colData(spe) containing the cell type identifiers.
#' Set to cell_type by default.
#' @param k Integer value. The number of distinct color to be generated, default is 30.
#' @param ngrid Integer value. The threshold (number of grid) used to find large ROIs.
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
#' coi <- c("Breast cancer", "Fibroblasts")
#' 
#' spe <- gridDensity(spe, coi = coi)
#' 
#' spe <- findROI(spe, coi = coi, method = "walktrap", steps = 5)
#' 
#' plotROI(spe, ngrid = 30, size = 0.3, alpha = 0.2)
#' 
plotROI <- function(spe, 
                    id = "cell_type", k = 30, 
                    ngrid = 1, showlegend = FALSE, ...){

  rois <- spe@metadata$roi

  coi <- spe@metadata$coi
  coi_clean <- janitor::make_clean_names(coi)

  dat <- as.data.frame(colData(spe))
  
  if (!is.null(coi)) {
    spe <- spe[, dat[, id] %in% coi]
  }
  
  posdat <- as.data.frame(spatialCoords(spe))
  
  dat <- as.data.frame(colData(spe)) |>
    cbind(posdat)
  
  set.seed(100)
  col.p <- randomcoloR::distinctColorPalette(k)
  
  xlim <- c(min(posdat[,1]), max(posdat[,1]))
  ylim <- c(min(posdat[,2]), max(posdat[,2]))
  plot.xlim <- xlim + c(-1e-10, 1e-10)
  plot.ylim <- ylim + c(-1e-10, 1e-10)
  
  filtered <- which(table(rois$component) >= ngrid)
  rois_filtered <- rois[rois$component %in% filtered, ]
  
  roi_plot <- plotSpatial(spe, ...) +
    geom_tile(data = rois_filtered, aes(x = xcoord, y = ycoord, fill = component), alpha = 0.6) +
    #scale_fill_manual(values = col.p) +
    scale_fill_manual(values = rep(col.p[1:10], 100)) +
    scale_x_continuous(limits = plot.xlim) +
    scale_y_continuous(limits = plot.ylim)
  
  if (isFALSE(showlegend)){
    roi_plot <- roi_plot +
      theme(legend.position = "none")
  } 

  return(roi_plot)
}

