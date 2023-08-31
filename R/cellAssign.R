#' Annotate all cells with contour level of cell type-specific density.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' @param name_to Colname in colData(spe) to store the annotation. 
#' @param NA_level Label for cells not falling in any of the regions. Default to 0.
#' @param levels Factor levels. 
#'
#' @return A SpatialExperiment object. An extra column is added to the colData.
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
#' spe <- contourAssign(spe, coi = coi)
#'
cellAssign <- function(spe, coi = NULL, 
                       name_to = paste0("at_level_",
                                           janitor::make_clean_names(coi)),
                          NA_level = "0", levels = NULL,
                          assign = c("contour", "roi"), ngrid = 20){
  
  if(length(assign) == 2L){
    assign <- "contour"
  } else if(!(assign %in% c("contour", "roi"))){
    stop("Parameter assign can only be either contour or roi.")
  }
  
  if(assign == "contour"){
    if(is.null(coi)){
      stop("Parameter coi is required for assigning cells to contour levels.")
    }
    all_areas <- getContourRegions(spe, coi = coi)
  } else {
    name_to <- "at_rois"
    NA_level = "no_roi"
    all_areas <- grid2sf(spe, ngrid = ngrid)
  }
  
  spe <- cellsInRegion(spe, all_areas, name_to = name_to,
                       NA_level = NA_level, levels = levels)
  
  return(spe)
}