#' Combine grids in each ROI to a sf region
#'
#' @param spe A SpatialExperiment object.
#' @param ngrid Integer. The threshold (minimum number of grids) used to filter small ROIs. Default to 20.
#'
#' @return List of ROIs saved as sf objects. 
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
#' spe <- findROI(spe, coi = coi, probs = 0.85)
#' 
#' rois <- grid2sf(spe, ngrid = 20)
#' 

grid2sf <- function(spe, ngrid = 20) {
  
  if (is.null(spe@metadata$roi))
    stop("ROI not yet computed!")
  
  rois <- as.data.frame(spe@metadata$roi)
  filtered <- which(table(rois$component) >= ngrid)
  rois_filtered <- rois[rois$component %in% filtered, ]
  
  grid_width <- spe@metadata$grid_info$xstep
  grid_height <- spe@metadata$grid_info$ystep
  canvas <- data.frame(x = rep(spe@metadata$grid_info$xlim, each = 2), 
                       y = rep(spe@metadata$grid_info$ylim, 2)) |>
    st_as_sf(coords = c("x", "y"))
  grids <- st_make_grid(canvas, n = c(spe@metadata$grid_info$dims[1], 
                                      spe@metadata$grid_info$dims[2]), what = "polygons")
  
  rois_sf <- list()
  for (rr in unique(rois_filtered$component)) {
    this_roi <- rois_filtered[rois_filtered$component == rr, ]
    this_roi_centroids <- st_as_sf(this_roi, coords = c("xcoord", "ycoord"))
    kp <- st_intersects(this_roi_centroids, grids)
    rois_sf[[rr]] <- st_as_sf(st_union(grids[as.numeric(kp)]))
  }
  
  return(rois_sf)
  
}