#' Plot ROIs on spatial.
#'
#' @param spe A SpatialExperiment object.
#' @param id Character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default.
#' @param show.legend Logical. Show legend or not.
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
#' plotROI(spe, size = 0.3, alpha = 0.2)
#'
plotROI <- function(spe,
                    id = "cell_type",
                    show.legend = FALSE, ...) {
  if (is.null(spe@metadata$roi)) {
    stop("ROI not yet computed!")
  }

  rois <- as.data.frame(spe@metadata$roi)

  coi <- spe@metadata$coi
  coi_clean <- janitor::make_clean_names(coi)

  dat <- as.data.frame(spe@colData)

  if (!is.null(coi)) {
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
  sf <- grid2sf(spe)
  rois_center <- do.call(rbind, lapply(sf, function(rr) {
    center <- sf::st_point_on_surface(rr)
    as.data.frame(sf::st_coordinates(center))
  }))

  rois_center <- as.data.frame(rois_center) |>
    rownames2col("component")

  roi_plot <- plotSpatial(spe, ...) +
    geom_tile(
      data = rois, aes(x = xcoord, y = ycoord, fill = component),
      alpha = 0.6
    ) +
    # scale_fill_manual(values = col.p) +
    annotate("text",
      x = rois_center$X, y = rois_center$Y,
      label = rois_center$component, color = "black", fontface = 2
    ) +
    scale_fill_manual(values = col.p) +
    scale_x_continuous(limits = plot.xlim) +
    scale_y_continuous(limits = plot.ylim)

  if (isFALSE(show.legend)) {
    roi_plot <- roi_plot +
      theme(legend.position = "none")
  }

  return(roi_plot)
}

utils::globalVariables(c("xcoord", "ycoord", "component"))
