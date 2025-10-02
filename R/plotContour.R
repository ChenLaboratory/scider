#' Plot contour lines.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs). 
#' All cell types are chosen if NULL or 'overall'.
#' @param overlay Character vector. Options are 'cell' (plot overlay on cells),
#' 'density' (overlay on density), or 'none'. Default to 'cell'.
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to 'cell_type' by default.
#' @param sub.level Character vector. Subset on specific level.
#' @param line.type shape of contour. See 'ggplot2::geom_path()'.
#' @param line.width size of contour.
#' @param line.alpha alpha of contour between 0 and 1.
#' @param reverseY Logical. Whether to reverse Y coordinates. Default is TRUE 
#' if the spe contains an image (even if not plotted) and FALSE if otherwise.
#' @param ... Aesthetic mappings to pass to \link[scider]{plotSpatial}, 
#' \link[scider]{plotDensity}, or \link[scider]{plotImage}, depending on 
#' the overlay.
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
#' plotContour(spe, coi = coi, line.width = 0.3, pt.alpha = 0.2)
#'
plotContour <- function(spe,
                        coi = NULL,
                        overlay = c("cell", "density", "none"),
                        id = "cell_type",
                        sub.level = NULL, 
                        line.type = 1,
                        line.width = 0.5,
                        line.alpha = 1,
                        reverseY = NULL,
                        ...) {
  coi_clean <- `if`(is.null(coi),"overall",cleanName(coi))
  contour_name <- paste(c(coi_clean,"contour"), collapse="_")

  if (!contour_name %in% names(spe@metadata)) {
    stop("Contour of interest doesn't exist. Please run getContour() first!")
  }
  
  contour_data <- as.data.frame(spe@metadata[[contour_name]])
  levs_drawn <- sort(unique(contour_data$level))
  levs_legend <- if (0 %in% levs_drawn) levs_drawn else c(0, levs_drawn)
  
  overlay <- overlay[1]
  if (overlay == "cell") {
    sub <- TRUE
    if (contour_name != "overall_roi" &&
        id %in% colnames(spe@colData) &&
        any(spe@colData[,id] %in% coi)) {
      sub <- spe@colData[[id]] %in% coi
    }
    p <- plotSpatial(spe[, sub],reverseY=reverseY,...)
  } else if (overlay == "density") {
    p <- plotDensity(spe, coi = coi,reverseY=reverseY,...)
  } else if (overlay == "none") {
    p <- plotImage(spe,reverseY=reverseY,...)
  } else {
    stop("Invalid 'overlay'.")
  }
  col.p <- grDevices::colorRampPalette(col.spec)(length(levs_legend))
  col.p <- rev(col.p)
  names(col.p) <- levs_legend
  
  if (is.null(sub.level)) {
    suppressMessages(p <- p +
                       ggplot2::geom_path(
                         data = contour_data,
                         ggplot2::aes(
                           x = x, y = y, group = group,
                           color = factor(level,
                                          levels = levs_legend)
                         ),
                         linewidth=line.width,
                         linetype=line.type,
                         alpha=line.alpha
                       ) +
                       ggplot2::scale_color_manual(name = "Density level", values = col.p,
                                                   breaks = levs_legend, drop = FALSE))
  } else {
    if (length(sub.level) == 1L & sub.level %in% contour_data$level) {
      suppressMessages(p <- p +
                         ggplot2::geom_path(
                           data = contour_data,
                           ggplot2::aes(
                             x = x, y = y, group = group,
                             color = level == sub.level
                           ),
                           linewidth=line.width,
                           linetype=line.type,
                           alpha=line.alpha
                         ) +
                         ggplot2::scale_color_manual(
                           name = paste0("level", sub.level, " density"),
                           values = c("royalblue", "tomato2")
                         ))
    } else {
      stop("The length sub.level is expected to be 1 and
           should be included in contour_data$level.")
    }
  }
  p <- update_bound(p,
                    x=SpatialExperiment::spatialCoords(spe)[, 1],
                    y=SpatialExperiment::spatialCoords(spe)[, 2])
  p <- p +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "x", y = "y") +
    ggplot2::ggtitle(paste(coi_clean, collapse=", "))
  return(p)
}

utils::globalVariables(c("x", "y", "group", "level"))
