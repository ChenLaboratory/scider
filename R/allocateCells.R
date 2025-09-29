#' Annotate all cells with contour level of cell type-specific density.
#'
#' @param spe A SpatialExperiment object.
#' @param to.roi Logical. Whether to allocate cells to ROIs.
#' @param roi Character. The name of the group or cell type on which
#' the roi is computed. If NULL, then the cell allocation will 
#' be performed for all detected roi Default to NULL.
#' @param to.contour Logical. Whether to allocate cells to contour levels.
#' @param contour Character. The name of the group or cell type on which
#' the contour level is computed. If NULL, then the cell allocation will 
#' be performed for all detected contours. Default to NULL.
#' 
#' @return A SpatialExperiment object. An extra column is added to the colData.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- gridDensity(spe)
#' coi <- "Breast cancer"
#' spe <- findROI(spe, coi = coi)
#' spe <- getContour(spe, coi = coi)
#' spe <- allocateCells(spe, contour = coi)
#'
allocateCells <- function(spe,
                          to.roi = TRUE,
                          roi = NULL,
                          to.contour = TRUE,
                          contour = NULL) {
    if (to.roi) {
        if (is.null(roi)) {
            roi_clean <- grep("_roi$", names(spe@metadata),value = TRUE)
        } else {
            roi_clean <- cleanName(roi)
            roi_clean <- paste(c(roi_clean,"roi"), collapse="_")
            if (!roi_clean %in% names(spe@metadata)) roi_clean <- NULL
        }
        
        if(length(roi_clean) == 0){
            message("No roi detected. Proceed without roi.")
        } else {
            for (r in roi_clean) {
                message(paste(
                    "Assigning cells to ROIs defined by",
                    janitor::make_clean_names(r,case="sentence",replace = c("_roi$"="")), 
                    "\n"
                ))
                rois <- spe@metadata[[r]]
                sf <- grid2sf(spe, rois$x,rois$y)
                # Unioning sf polygons with same ROIs
                all_areas <- lapply(unique(rois$component), function(xx) {
                    sf::st_as_sf(sf::st_union(sf::st_sfc(sf[rois$component == xx]),is_coverage = TRUE))
                })
                names(all_areas) <- unique(rois$component)
                
                spe <- cellsInRegion(spe, all_areas,
                                     name_to = r,
                                     NA_level = "None", levels = NULL
                )
            }
        }
    }
    
    if (to.contour) {
        if (is.null(contour)) {
            contour_clean <- grep("_contour$", names(spe@metadata),value = TRUE)
        } else {
            contour_clean <- cleanName(contour)
            contour_clean <- paste(c(contour_clean,"contour"), collapse="_")
            if (!contour_clean %in% names(spe@metadata)) contour_clean <- NULL
        }
      
        if (length(contour_clean)==0) {
            message("No contour detected. Proceed without contour.")
        } else {
            for (cc in contour_clean) {
                message(paste(
                    "Assigning cells to contour levels of",
                    janitor::make_clean_names(cc,case="sentence",replace = c("_contour$"="")), 
                    "\n"
                ))
            
            all_areas <- getContourRegions(spe, contour_name = cc)
            
            spe <- cellsInRegion(spe, all_areas, name_to = cc,
                                 NA_level = 0, levels = NULL)
            }
        }
    }
    
    return(spe)
}
