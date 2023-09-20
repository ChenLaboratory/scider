#' Combine grids in each ROI to a sf region
#'
#' @param spe A SpatialExperiment object.
#'
#' @return List of ROIs saved as sf objects.
#'
grid2sf <- function(spe) {
    if (is.null(spe@metadata$roi)) {
        stop("ROI not yet computed!")
    }

    rois <- as.data.frame(spe@metadata$roi)

    grid_width <- spe@metadata$grid_info$xstep
    grid_height <- spe@metadata$grid_info$ystep
    canvas <- data.frame(
        x = rep(spe@metadata$grid_info$xlim, each = 2),
        y = rep(spe@metadata$grid_info$ylim, 2)
    ) |>
        sf::st_as_sf(coords = c("x", "y"))
    grids <- sf::st_make_grid(canvas, n = c(
        spe@metadata$grid_info$dims[1],
        spe@metadata$grid_info$dims[2]
    ), what = "polygons")

    rois_sf <- list()
    for (rr in unique(rois$component)) {
        this_roi <- rois[rois$component == rr, ]
        this_roi_centroids <- sf::st_as_sf(this_roi,
            coords = c("xcoord", "ycoord")
        )
        kp <- sf::st_intersects(this_roi_centroids, grids)
        rois_sf[[rr]] <- sf::st_as_sf(sf::st_union(grids[as.numeric(kp)]))
    }

    return(rois_sf)
}
