#' Check which cells are in which regions
#'
#' @param spe A SpatialExperiment object.
#' @param sf List or an sf object that represents a region or an ROI. 
#' @param name_to Colname in colData(spe) to store the annotation. 
#' @param NA_level Label for cells not falling in any of the regions. 
#' @param levels Factor levels. 
#'
#' @return An spe object. 
#' @export
#'
#' @examples
cellsInRegion <- function(spe, sf, name_to, NA_level = NULL, levels = NULL) {
  
  if (length(sf) > 1L) {
    sf_classes <- sapply(sf, class)[1, ]
  }
  if (length(sf) == 1L) {
    sf_classes <- class(sf)[1]
  }
  if (any(sf_classes != "sf"))
    stop("One or more regions not converted to the sf class!")
  
  # all cells
  xy_allcells <- st_as_sf(as.data.frame(spatialCoords(spe)), coords = c("x_centroid", "y_centroid"))
  
  # calculate overlaps
  isIn <- list()
  for (aa in names(sf)) {
    # contour region
    this_area <- sf[[aa]]
    # calculate intersection
    overlap_ind <- st_intersects(xy_allcells, this_area, sparse = FALSE)
    overlap_ind <- which(overlap_ind == 1)
    isIn[[aa]] <- overlap_ind
  }
  
  # annotate colData
  to_append <- rep(NA_character_, nrow(colData(spe)))
  for (aa in names(isIn)) {
    to_append[isIn[[aa]]] <- aa
  }
  if (anyNA(to_append)) {
    if (is.null(NA_level))
      stop("Need to specify `NA_level` as labels for cells not in any of the regions!")
    to_append[is.na(to_append)] <- NA_level
  }
  
  if (is.null(levels))
    levels <- names(table(to_append))
  colData(spe)[[name_to]] <- factor(to_append, levels = levels)
  
  return(spe)
  
}