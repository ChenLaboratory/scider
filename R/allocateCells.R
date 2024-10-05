#' Annotate all cells with contour level of cell type-specific density.
#'
#' @param spe A SpatialExperiment object.
#' @param to.roi Logical. Whether to allocate cells to ROIs.
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
#' spe <- getContour(spe, coi=coi)
#' spe <- allocateCells(spe, contour)
#'
allocateCells <- function(
    spe,
    to.roi = TRUE,
    to.contour = TRUE,
    contour = NULL) {
    if (to.roi) {
        if (is.null(spe@metadata$roi)) {
            message("No ROI detected.")
        } else {
            message(paste(
                "Assigning cells to ROIs defined by",
                paste(spe@metadata$coi, collapse = ", "), "\n"
            ))
            all_areas <- grid2sf(spe)
            name_to <- "roi"
            NA_level <- "None"
            spe <- cellsInRegion(spe, all_areas,
                name_to = name_to,
                NA_level = NA_level, levels = NULL
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
                if(length(contour_clean)>1)
                    contour_clean <- paste(sort(contour_clean), collapse="_")
                if(! paste0(contour_clean, "_contour") %in% names(spe@metadata)){
                    message("Specified contour not detected. Proceed without contour.")
                    ind <- integer(0)
                } else {
                    ind <- grep(paste0(contour_clean, "_contour"), names(spe@metadata))
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
                    NA_level <- 0

                    spe <- cellsInRegion(spe, all_areas, name_to = name_to,
                        NA_level = NA_level, levels = NULL)
                }
            }
        }
    }

    return(spe)
}
