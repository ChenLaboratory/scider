#' Check which cells are in which regions
#'
#' @param spe A SpatialExperiment object.
#' @param region List or an sf object that represents a region or an ROI.
#' @param name_to Colname in colData(spe) to store the annotation.
#' @param NA_level Label for cells not falling in any of the regions.
#' Default to 0.
#' @param levels Factor levels.
#'
#' @return A SpatialExperiment object. The region information of each cell is
#' stored in the colData.
#'
#'
cellsInRegion <- function(spe, region, name_to,
                          NA_level = "0", levels = NULL) {
    if (is.null(names(region))) {
        warning("The region input is unnamed! We recommend a named list of
            region object(s) as input!")
    }
    if (is.null(NA_level)) {
        stop("Need to specify `NA_level` as labels for cells not
             in any of the regions!")
    }

    # all cells
    xy_allcells <- sf::st_as_sf(as.data.frame(SpatialExperiment::spatialCoords(spe)),
        coords = c("x_centroid", "y_centroid")
    )

    # calculate overlaps
    region_names <- names(region) %||% as.character(seq_along(region))
    region <- do.call(rbind,region)
    isIn <- sf::st_intersects(region,xy_allcells)

    # annotate colData
    to_append <- rep_len(NA_level, nrow(colData(spe)))
    for (i in seq_along(region_names)) {
        to_append[isIn[[i]]] <- region_names[i]
    }

    if (is.null(levels)) {
        # levels <- unique(to_append)[order(as.numeric(unique(to_append)))]
        val <- unique(to_append)[unique(to_append) != NA_level]
        levels <- c(NA_level, val[order(as.numeric(val))])
        colData(spe)[[name_to]] <- factor(to_append, levels = levels)
    }
    return(spe)
}
