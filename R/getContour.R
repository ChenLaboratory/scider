#' Get contour from density
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' @param bins An integer. Number of contour levels.
#' @param binwidth A numeric scale of the smoothing bandwidth.
#' @param breaks A numeric scale referring to the breaks in
#' `ggplot2:::contour_breaks`.
#'
#' @return A SpatialExperiment object. An sf object of the contour region of
#' the specified level is stored in the metadata of the
#' SpatialExperiment object.
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
getContour <- function(spe, coi, bins = NULL,
                       binwidth = NULL, breaks = NULL) {
    if (is.null(spe@metadata$grid_density)) {
        stop("Have to calculate grid density, run gridDensity() first!")
    }

    dens <- spe@metadata$grid_density
    dups <- duplicated(dens[, c("y_grid", "x_grid"), drop = FALSE],
        fromLast = TRUE
    )
    dens <- dens[!dups, , drop = FALSE]

    coi_clean <- janitor::make_clean_names(coi)
    dens_cols <- paste("density", coi_clean, sep = "_")

    if (!all(dens_cols %in% colnames(dens))) {
        stop("Density of COI is not yet computed.")
    }

    if (length(dens_cols) > 1L) {
        message("Plotting contour of average density of input COIs. ")
        dens$density_coi_average <- rowMeans(dens[, which(colnames(dens) %in%
            dens_cols),
        drop = FALSE
        ])
    } else {
        dens$density_coi_average <- dens[, dens_cols]
    }
    dens <- dens[, c(seq_len(5), which(colnames(dens) ==
        "density_coi_average"))]

    # levels for contour
    if (is.null(bins) && is.null(binwidth) && is.null(breaks)) {
        message("Using bins = 10 to draw contours.")
        bins <- 10L
    }
    if (!is.null(bins)) binwidth <- breaks <- NULL
    if (is.null(bins) && !is.null(binwidth)) breaks <- NULL

    # filter out negative densities when calculating contours
    dens <- dens[dens$density_coi_average > 0L, ]

    # note that when calculating contours, density is not filtered at any
    # quantile cutoff!
    contour <- compute_group(dens,
        z.range = range(dens$density_coi_average, na.rm = TRUE, finite = TRUE),
        bins = bins,
        binwidth = binwidth,
        breaks = breaks,
        na.rm = FALSE
    )

    contour$level <- as.factor(as.numeric(as.factor(contour$cutoff)))
    spe@metadata[[paste(coi_clean,
        "contour",
        sep = "_"
    )]] <- S4Vectors::DataFrame(contour)

    return(spe)
}



#### lower level functions for computing the contours.
#### NEED TO CLEAN UP LATER!!!
# compute contour groups (grabbed from ggplot2)
xyz_to_isolines <- function(data, breaks) {
    isoband::isolines(
        x = sort(unique00(data$x_grid)),
        y = sort(unique00(data$y_grid)),
        z = isoband_z_matrix(data),
        levels = breaks
    )
}

isoband_z_matrix <- function(data) {
    # Convert vector of data to raster
    x_pos <- as.integer(factor(data$x_grid,
        levels =
            sort(unique00(data$x_grid))
    ))
    y_pos <- as.integer(factor(data$y_grid,
        levels =
            sort(unique00(data$y_grid))
    ))
    nrow <- max(y_pos)
    ncol <- max(x_pos)
    raster <- matrix(NA_real_, nrow = nrow, ncol = ncol)
    raster[cbind(y_pos, x_pos)] <- data$density_coi_average
    raster
}

iso_to_path <- function(iso, group = 1) {
    lengths <- vapply(iso, function(x) length(x$x), integer(1))

    if (all(lengths == 0)) {
        message("Zero contours were generated.")
        return(data_frame00())
    }

    levels <- names(iso)
    xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
    ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
    ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
    item_id <- rep(seq_along(iso), lengths)

    # Add leading zeros so that groups can be properly sorted
    groups <- paste(group, sprintf("%03d", item_id), sprintf("%03d", ids),
        sep = "-"
    )
    groups <- factor(groups)

    data_frame00(
        level = rep(levels, lengths),
        x = xs,
        y = ys,
        piece = as.integer(groups),
        group = groups,
        .size = length(xs)
    )
}

compute_group <- function(data, z.range, bins = NULL,
                          binwidth = NULL, breaks = NULL, na.rm = FALSE) {
    breaks <- contour_brks(z.range, bins, binwidth, breaks)
    isolines <- xyz_to_isolines(data, breaks)
    path_df <- iso_to_path(isolines, data$group[1])
    path_df$cutoff <- as.numeric(path_df$level)
    # path_df$level <- as.numeric(path_df$level)
    # path_df$nlevel <- scales::rescale_max(path_df$level)
    path_df
}
