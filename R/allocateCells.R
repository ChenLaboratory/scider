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
            roi_clean <- gsub("_roi$", "", roi)
            roi_clean <- janitor::make_clean_names(roi_clean)
            roi_clean <- paste(c(sort(roi_clean),"roi"), collapse="_")
        }
        
        if(length(roi_clean) == 0 || !roi_clean %in% names(spe@metadata)){
            message("No roi detected. Proceed without roi.")
            roi_clean <- integer(0)
        }
        
        for (r in roi_clean) {
            message(paste(
                "Assigning cells to ROIs defined by",
                janitor::make_clean_names(r,case="sentence",replace = c("roi"="")), 
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
    
    if (to.contour) {
        ind <- grep("_contour", names(spe@metadata))
        coi_2 <- NULL
        if (length(ind) == 0) {
            message("No contour detected.")
        } else {
            if (!is.null(contour)){
                contour_clean <- janitor::make_clean_names(contour)
                contour_clean <- paste(c(sort(contour_clean),"contour"), collapse="_")
                if(! contour_clean %in% names(spe@metadata)){
                    message("Specified contour not detected. Proceed without contour.")
                    ind <- integer(0)
                } else {
                    ind <- grep(contour_clean, names(spe@metadata))
                    coi_2 <- paste(contour, collapse=", ")
                }
            }
            
            for (i in ind) {
                coi <- janitor::make_clean_names(names(spe@metadata)[i],
                                                 case = "sentence", replace = c("contour" = ""))
                if(!is.null(coi_2)) coi <- contour
                
                if(all(paste0("density_", janitor::make_clean_names(coi)) %in% colnames(spe@metadata$grid_density))){
                    message(paste(
                        "Assigning cells to contour levels of",
                        paste(coi, collapse=", "), "\n"
                    ))
                    all_areas <- getContourRegions(spe, coi = coi)
                    name_to <- names(spe@metadata)[i]
                    
                    spe <- cellsInRegion(spe, all_areas, name_to = name_to,
                                         NA_level = 0, levels = NULL)
                }
            }
        }
    }
    
    return(spe)
}
