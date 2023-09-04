#' Check which cells are in which regions
#'
#' @param spe A SpatialExperiment object.
#' @param region List or an sf object that represents a region or an ROI. 
#' @param name_to Colname in colData(spe) to store the annotation. 
#' @param NA_level Label for cells not falling in any of the regions. Default to 0.
#' @param levels Factor levels. 
#'
#' @return A SpatialExperiment object. The region information of each cell is stored in the colData. 
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
#' region <- getContourRegions(spe, coi = coi)
#' 
#' spe <- cellsInRegion(spe, region, name_to = "breast_cancer_contour_level")
#' 

cellsInRegion <- function(spe, region, name_to, NA_level = "0", levels = NULL) {
  
  if (length(region) > 1L) {
    sf_classes <- sapply(region, class)[1, ]
  }
  if (length(region) == 1L) {
    sf_classes <- class(region)[1]
  }
  if (any(sf_classes != "sf"))
    stop("One or more regions not converted to the sf class!")
  
  if (is.null(names(region)))
    warning("The region input is unnamed! We recommend a named list of region object(s) as input!")
    
  # all cells
  xy_allcells <- sf::st_as_sf(as.data.frame(spatialCoords(spe)), coords = c("x_centroid", "y_centroid"))
  
  # calculate overlaps
  isIn <- list()
  for (aa in 1:length(region)) {
    # contour region
    this_area <- region[[aa]]
    # calculate intersection
    overlap_ind <- sf::st_intersects(xy_allcells, this_area, sparse = FALSE)
    overlap_ind <- which(overlap_ind == 1)
    isIn[[aa]] <- overlap_ind
  }
  
  if (!is.null(names(region))) {
    names(isIn) <- names(region)
  } else {
    names(isIn) <- as.character(1:length(region))
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
    levels <- unique(to_append)
  colData(spe)[[name_to]] <- factor(to_append, levels = levels)
  
  return(spe)
  
}