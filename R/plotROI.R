#' Plot ROIs on spatial.
#'
#' @param spe A SpatialExperiment object.
#' @param roi Character. The name of the group or cell type on which
#' the roi is computed. All cell types are chosen if NULL or 'overall'.
#' @param id Character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default.
#' @param label Logical. Show ROI label or not.
#' @param show.legend Logical. Show legend or not.
#' @param ... Parameters pass to \link[scider]{plotSpatial}
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
#' plotROI(spe, roi = coi, pt.size = 0.3, pt.alpha = 0.2)
#'
plotROI <- function(spe,
                    roi = NULL,
                    id = "cell_type",
                    label = TRUE,
                    show.legend = FALSE, 
                    ...) {
  roi_clean <- `if`(is.null(roi),"overall",cleanName(roi))
  roi_clean_name <- paste(c(roi_clean,"roi"), collapse="_")
  
  if (is.null(spe@metadata[[roi_clean_name]])) {
    stop("ROI of interest doesn't exist. Please run findROI(), or specify 'roi'.")
  }
  
  rois <- as.data.frame(spe@metadata[[roi_clean_name]])
  
  # Filter background points to roi only.
  if (roi_clean_name != "overall_roi" &&
      id %in% colnames(spe@colData) &&
      any(spe@colData[,id] %in% roi)) {
    spe <- spe[, spe@colData[,id] %in% roi]
  }
  
  nROIs <- nlevels(rois$component)
  col.p <- selectColor(nROIs)
  
  # Get centres of ROIS
  sf <- grid2sf(spe, rois$x,rois$y)
  sf <- lapply(unique(rois$component), function(xx) {
    sf::st_union(sf::st_sfc(sf[rois$component == xx]))
  })
  names(sf) <- as.character(unique(rois$component))
  rois_center <- do.call(rbind, lapply(sf, function(rr) {
    center <- sf::st_point_on_surface(rr)
    as.data.frame(sf::st_coordinates(center))
  }))
  rois_center <- as.data.frame(rois_center)
  rois_center$component <- names(sf)
  
  # Turn sf into a dataframe for geom_polygon
  # Dealing with MULTIPOLYGON (which can happen with mergeROI)
  sf <- lapply(sf, function(xx) sf::st_cast(xx,"POLYGON"))
  poly_groups <- list(rep.int(1:length(sf),lengths(sf)))
  sf <- as.data.frame(sf::st_coordinates(do.call(c,sf)))
  sf$component <- rep.int(rois_center$component,
                          stats::aggregate(
                            tabulate(sf$L2),
                            by=poly_groups,
                            FUN=sum)[[2]])
  # Plotting
  p <- plotSpatial(spe,...) +
    geom_polygon(
      data = sf,
      aes(
        x=X,
        y=Y,
        group = L2,
        fill = component,
        alpha = 0.6,
        subgroup = L1
      ),color=NA,
      inherit.aes = F) +
    scale_fill_manual(values = col.p) +
    ggtitle(paste0("ROI (", paste(roi_clean, collapse=", "), ")"))
  p <- update_bound(p, y=sf$Y)

  if (isFALSE(show.legend)) {
    p <- p +
      theme(legend.position = "none")
  }
  
  if (label) {
    p <- p +
      annotate("text",
               x = rois_center$X, y = rois_center$Y,
               label = rois_center$component, color = "black", fontface = 2
      )
  }
  
  return(p)
}


utils::globalVariables(c("xcoord", "ycoord", "component"))
