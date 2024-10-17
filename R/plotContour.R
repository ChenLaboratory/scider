#' Plot contour lines.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs). 
#' All cell types are chosen if NULL or 'overall'.
#' @param overlay Character vector. Options are 'cell' (plot overlay on cells),
#' 'density' (overlay on density), or 'none'. Default to 'cell'.
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to 'cell_type' by default.
#' @param sub.level Character vector. Subset on specific level.
#' @param line.type shape of contour. See 'ggplot2::geom_path()'.
#' @param line.width size of contour.
#' @param line.alpha alpha of contour between 0 and 1.
#' @param ... Aesthetic mappings to pass to 'plotSpatial()' or 
#' 'plotDensity()', depending on the overlay.
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
#' coi <- "Breast cancer"
#'
#' spe <- getContour(spe, coi = coi)
#'
#' plotContour(spe, coi = coi, line.width = 0.3, pt.alpha = 0.2)
#'
plotContour <- function(spe,
                        coi = NULL,
                        overlay = c("cell", "density", "none"),
                        id = "cell_type",
                        sub.level = NULL, 
                        line.type = 1,
                        line.width = 0.5,
                        line.alpha = 1,
                        
                        ...) {
    if ( !is.null(coi) & !("overall" %in% coi) ){
        if ( ! all(coi %in% names(table(colData(spe)[[id]]))) ) {
            stop("coi not in colData(spe)[[id]]!")
        }
    } else coi <- "overall"

    coi_clean <- janitor::make_clean_names(coi)
    if (length(coi_clean) > 1L) coi_clean <- paste(sort(coi_clean), collapse="_")
    coi_clean_contour <- paste(coi_clean, "contour", sep = "_")

    if (!coi_clean_contour %in% names(spe@metadata)) {
        stop("Contour of interest doesn't exist. Please run getContour() first!")
    }

    contour_data <- as.data.frame(spe@metadata[[coi_clean_contour]])
    levs <- unique(contour_data$level)
    nlevs <- length(levs)

    overlay <- overlay[1]
    if (overlay == "cell") {
        sub <- TRUE
        if(all(coi != "overall"))
            sub <- colData(spe)[[id]] %in% coi
        p <- plotSpatial(spe[, sub], ...)
    } else if (overlay == "density") {
        p <- plotDensity(spe, coi = coi, ...)
    } else if (overlay == "none") {
        p <- plotSpatial(spe[, FALSE], ...)
    } else {
        stop("Invalid 'overlay'.")
    }

    col.p <- grDevices::colorRampPalette(col.spec)(
        length(unique(contour_data$level)))

    if (is.null(sub.level)) {
        suppressMessages(p <- p +
            ggplot2::geom_path(
                data = contour_data,
                ggplot2::aes(
                    x = x, y = y, group = group,
                    color = level
                ),
                linewidth=line.width,
                linetype=line.type,
                alpha=line.alpha
            ) +
            scale_color_manual(name = "Density level", values = rev(col.p)))
    } else {
        if (length(sub.level) == 1L & sub.level %in% contour_data$level) {
            suppressMessages(p <- p +
                ggplot2::geom_path(
                    data = contour_data,
                    ggplot2::aes(
                        x = x, y = y, group = group,
                        color = level == sub.level
                    ),
                    linewidth=line.width,
                    linetype=line.type,
                    alpha=line.alpha
                ) +
                scale_color_manual(
                    name = paste0("level", sub.level, " density"),
                    values = c("royalblue", "tomato2")
                ))
        } else {
            stop("The length sub.level is expected to be 1 and
           should be included in contour_data$level.")
        }
    }

    p <- p +
        theme_classic() +
        labs(x = "x", y = "y") +
        ggtitle(paste(coi, collapse=", "))
    return(p)
}


utils::globalVariables(c("x", "y", "group", "level"))
