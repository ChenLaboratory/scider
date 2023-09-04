#' Annotate all cells with contour level of cell type-specific density.
#'
#' @param spe A SpatialExperiment object.
#' @param to.roi Logical. Whether to allocate cells to ROIs.
#' @param to.contour Logical. Whether to allocate cells to contour levels.
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
#' spe <- findROI(spe, coi = coi)
#' 
#' spe <- allocateCells(spe)
#'

allocateCells <- function(spe, 
  to.roi = TRUE,
  to.contour = TRUE){

  if(to.roi){
    if(is.null(spe@metadata$roi)){ 
      message("No ROI detected.")
    } else {
      cat(paste("Assigning cells to ROIs defined by", paste(spe@metadata$coi, collapse=", "), "\n"))
      all_areas <- grid2sf(spe)
      name_to <- "roi"
      NA_level <- "None"
      spe <- cellsInRegion(spe, all_areas, name_to = name_to,
                       NA_level = NA_level, levels = NULL)
    }
  }
  
  if(to.contour){
    ind <- grep("_contour", names(spe@metadata))
    if(length(ind) == 0){
      message("No contour detected.")
    } else {
      for(i in ind){
        coi <- janitor::make_clean_names(names(spe@metadata)[i], 
                case = "sentence", replace = c("contour" = ""))
        cat(paste("Assigning cells to contour levels of", coi, "\n"))
      
        all_areas <- getContourRegions(spe, coi = coi)

        name_to <- names(spe@metadata)[i]
        NA_level <- 0

        spe <- cellsInRegion(spe, all_areas, name_to = name_to,
                       NA_level = NA_level, levels = NULL)
      }
    }
  }

  return(spe)
}
