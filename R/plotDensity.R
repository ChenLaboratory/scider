#' Plot grid-based density.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs) to be plotted.
#' Default to all cell types.
#' @param probs Numeric value between 0 and 1, used for filtering
#' uninformative grid, default is 0.5.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe <- gridDensity(spe)
#'
#' plotDensity(spe, coi = "Breast cancer")
#'
#' plotDensity(spe, coi = "Fibroblasts")
#'
plotDensity <- function(spe, coi = NULL, probs = 0.5) {
    grid_data <- as.data.frame(spe@metadata$grid_density)

    if (is.null(coi)) coi <- "overall"
    if (length(coi) >= 2) coi <- coi[coi!="overall"]
    coi_clean <- janitor::make_clean_names(coi)

    dens_cols <- paste("density", coi_clean, sep = "_")

    if (!all(dens_cols %in% colnames(grid_data))) {
        stop("Density of COI is not yet computed.")
    }

    grid_data$density_coi_average <- rowSums(as.matrix(
        grid_data[, which(colnames(grid_data) %in% dens_cols),
            drop = FALSE]
    ))

    kp <- grid_data$density_coi_average >=
        quantile(grid_data$density_coi_average,
            probs = probs
        )

    xstep <- spe@metadata$grid_info$xstep
    ystep <- spe@metadata$grid_info$ystep

    p <- ggplot() +
        geom_tile(
            data = grid_data[kp, ],
            aes(
                x = x_grid, y = y_grid,
                fill = density_coi_average
            )
        ) + 
        coord_fixed() +
        theme_classic() +
        scale_fill_gradientn(colours = rev(col.spec)) +
        labs(x = "x", y = "y", fill = "Density") +
        lims(
            x = c(
                min(grid_data[, "x_grid"]) - xstep/2,
                max(grid_data[, "x_grid"]) + xstep/2
            ),
            y = c(
                min(grid_data[, "y_grid"]) - ystep/2,
                max(grid_data[, "y_grid"]) + ystep/2
            )
        ) +
    ggtitle(paste(coi, collapse=", "))

    return(p)
}

utils::globalVariables(c("x_grid", "y_grid", "density_coi_average"))
