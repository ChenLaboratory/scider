#' Get contour from density
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' @param equal.cell Logical. Whether to use produce contour levels so that
#' there are roughly the same number of cells of the COI at each level. 
#' @param bins An integer. Number of contour levels.
#' @param binwidth A numeric scale of the smoothing bandwidth.
#' @param breaks A numeric scale referring to the breaks in
#' `ggplot2:::contour_breaks`.
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default. Only needed when
#' \code{equal.cell = TRUE}. 
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
getContour <- function(spe, coi = NULL, equal.cell = FALSE, bins = NULL,
                       binwidth = NULL, breaks = NULL, id = "cell_type") {
    
    if (is.null(spe@metadata$grid_density)) {
        stop("Have to calculate grid density, run gridDensity() first!")
    }

    if (!id %in% colnames(colData(spe))) {
        stop(paste(id, "is not a column of the colData."))
    }

    if (is.null(coi)) {
        coi <- names(table(colData(spe)[[id]]))
    }

    if (length(which(!coi %in% names(table(colData(spe)[[id]])))) > 0L) {
        stop(paste(paste0(
            coi[which(!coi %in%
                names(table(colData(spe)[[id]])))],
            collapse = ", "
        ), "not found in data!", sep = " "))
    }

    coi_clean <- janitor::make_clean_names(coi)
    dens_cols <- paste("density", coi_clean, sep = "_")

    # grid level density data
    dens <- spe@metadata$grid_density
    dups <- duplicated(dens[, c("y_grid", "x_grid"), drop = FALSE],
        fromLast = TRUE
    )
    dens <- dens[!dups, , drop = FALSE]
    dens <- as.data.frame(dens)

    if (!all(dens_cols %in% colnames(dens))) {
        stop("Density of COI is not yet computed.")
    }

    if (length(dens_cols) > 1L) {
        message("Plotting contour of total density of input COIs. ")
        dens$density_coi <- rowSums(dens[, which(colnames(dens) %in%
            dens_cols),
        drop = FALSE
        ])
    } else {
        dens$density_coi <- dens[, dens_cols]
    }
    
    dens <- dens[, c(seq_len(5), which(colnames(dens) ==
        "density_coi"))]

    # filter out negative densities when calculating contours
    dens <- dens[dens$density_coi > 0L, ]

    # levels for contour
    if (!equal.cell) {
        if (is.null(bins) && is.null(binwidth) && is.null(breaks)) {
        message("Using bins = 10 to draw contours.")
        bins <- 10L
        }
        if (!is.null(bins)) binwidth <- breaks <- NULL
        if (is.null(bins) && !is.null(binwidth)) breaks <- NULL
    } else {
        # calculate density level breaks to get roughly the same number of cells at each level
        if (is.null(bins)) {
            message("Using bins = 10 to draw contours with equal cell numbers.")
            bins <- 10L
            binwidth <- breaks <- NULL
        } else {
            message(paste("Using bins =", bins, "to draw contours with equal cell numbers.", sep = " "))
        }
        ## count no of cells of coi in each grid
        ## note this no can be very different from the expected no
        coi_coords <- as.data.frame(spatialCoords(spe)[rownames(colData(spe))[colData(spe)[[id]] %in% coi], ])
        coi_coords$x_node <- vapply(coi_coords$x_centroid, function(xx) {
            which.min(abs(spe@metadata$grid_info$xcol - xx))
        }, numeric(1))
        coi_coords$y_node <- vapply(coi_coords$y_centroid, function(yy) {
            which.min(abs(spe@metadata$grid_info$yrow - yy))
        }, numeric(1))
        coi_coords$node <- paste(coi_coords$x_node, coi_coords$y_node, sep = "-")
        freq <- c(table(coi_coords$node))
        dens$n <- freq[dens$node]
        dens$n <- ifelse(is.na(dens$n), 0L, dens$n)
        dens_expanded <- rep(dens$density_coi, times = dens$n)
        dens_expanded <- dens_expanded[dens_expanded > 0L]
        qq <- seq(0, 1, round(1 / bins, 1))[-1]
        if (qq[length(qq)] == 1L) qq <- qq[-length(qq)]
        breaks <- as.vector(quantile(dens_expanded, probs = qq))
        binwidth <- bins <- NULL
    }

    # note that when calculating contours, density is not filtered at any
    # quantile cutoff!
    contour <- compute_group(dens,
        z.range = range(dens$density_coi, na.rm = TRUE, finite = TRUE),
        bins = bins,
        binwidth = binwidth,
        breaks = breaks,
        na.rm = FALSE
    )

    contour$level <- as.factor(as.numeric(as.factor(contour$cutoff)))

    coi_clean_output <- ifelse(length(coi_clean) == 1L, coi_clean, "coi")
    spe@metadata[[paste(coi_clean_output,
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
    raster[cbind(y_pos, x_pos)] <- data$density_coi
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
