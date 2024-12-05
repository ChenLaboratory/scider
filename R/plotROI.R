#' Plot ROIs on spatial.
#'
#' @param spe A SpatialExperiment object.
#' @param id Character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default.
#' @param label Logical. Show ROI label or not.
#' @param show.legend Logical. Show legend or not.
#' @param ... Aesthetic mappings pass for point.
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
#' plotROI(spe, pt.size = 0.3, pt.alpha = 0.2)
#'
plotROI <- function(spe,
                     id = "cell_type",
                     label = TRUE,
                     show.legend = FALSE, ...) {
  if (is.null(spe@metadata$roi)) {
    stop("ROI not yet computed!")
  }
  
  rois <- as.data.frame(spe@metadata$roi)
  
  coi <- spe@metadata$coi
  #coi_clean <- janitor::make_clean_names(coi)
  
  dat <- as.data.frame(spe@colData)
  
  if (!is.null(coi) && !("overall" %in% coi) && 
      (is.null(id) ||id %in% colnames(dat))) {
    spe <- spe[, dat[, id] %in% coi]
  }
  
  posdat <- as.data.frame(spatialCoords(spe))
  
  dat <- as.data.frame(spe@colData) |>
    cbind(posdat)
  
  nROIs <- nlevels(rois$component)
  col.p <- selectColor(nROIs)
  
  xlim <- spe@metadata$grid_info$xlim
  ylim <- spe@metadata$grid_info$ylim
  plot.xlim <- xlim + c(-1e-10, 1e-10)
  plot.ylim <- ylim + c(-1e-10, 1e-10)
  
  # filtered <- names(which(table(rois$component) >= ngrid))
  # rois_filtered <- as.data.frame(rois[rois$component %in% filtered, ])
  
  # for(n in colnames(colData(spe))){
  #  if (!(n %in% colnames(rois_filtered))){
  #    rois_filtered[, n] <- "dummy"
  #  }
  # }
  
  # Label ROI numbers at the center
  sf <- grid2sf(spe, rois$x,rois$y)
  sf = lapply(unique(rois$component), function(xx) {
    sf::st_union(sf::st_sfc(sf[rois$component == xx]))
  })
  
  rois_center <- do.call(rbind, lapply(sf, function(rr) {
    center <- sf::st_point_on_surface(rr)
    as.data.frame(sf::st_coordinates(center))
  }))
  
  rois_center <- as.data.frame(rois_center) |>
    rownames2col("component")
  
  # Plotting
  roi_plot = plotSpatial(spe, ...) + 
    geom_sf(
      data = sf::st_as_sfc(unlist(sf,recursive=F)),
      aes(
        fill = unique(rois$component),
        alpha = 0.6,
      ),color=NA,
      inherit.aes = F) +
    # scale_fill_manual(values = col.p) +
    scale_fill_manual(values = col.p) +
    scale_x_continuous(limits = plot.xlim) +
    scale_y_continuous(limits = plot.ylim) +
    ggtitle(paste0("ROI (", paste(coi, collapse=", "), ")"))
  
  if (isFALSE(show.legend)) {
    roi_plot <- roi_plot +
      theme(legend.position = "none")
  }

  if (label) {
    roi_plot <- roi_plot +
    annotate("text",
             x = rois_center$X, y = rois_center$Y,
             label = rois_center$component, color = "black", fontface = 2
    )
  }

  return(roi_plot)
}


utils::globalVariables(c("xcoord", "ycoord", "component"))
