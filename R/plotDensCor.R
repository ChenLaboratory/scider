#' Plot density correlation between two cell types
#'
#' @param spe A SpatialExperiment object.
#' @param celltype1 Cell type 1 to compare.
#' @param celltype2 Cell type 2 to compare.
#' @param by.roi Logical. Plot facet by ROIs or not. Default to TRUE.
#' @param probs A numeric scalar. The threshold of proportion that used to
#' filter grids by density when ROIs have not been identified previously.
#' Ignored if 'roi' is present in the 'metadata' component of spe. Default to 0.85.
#' @param fit Character. Options are "spline" and "linear".
#' @param df Integer. Degrees of freedom of the spline fit.
#' Default to 3 (i.e., a cubic spline fit).
#' @param ... aesthetic mappings to pass to `ggplot2::aes()`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' coi <- c("Breast cancer", "Fibroblasts")
#'
#' spe <- gridDensity(spe, coi = coi)
#'
#' spe <- findROI(spe, coi = coi, method = "walktrap")
#'
#' plotDensCor(spe, celltype1 = "Breast cancer", celltype2 = "Fibroblasts")
#'
plotDensCor <- function(spe, celltype1 = NULL, celltype2 = NULL,
                        by.roi = TRUE, probs = 0.85, 
                        fit = c("spline", "linear"), df = 3, ...) {
    if (!("grid_density" %in% names(spe@metadata))) {
        stop("Please run gridDensity before using this function.")
    }

    dens_dat <- as.data.frame(spe@metadata$grid_density)
    # clean names
    ct1 <- paste0("density_", janitor::make_clean_names(celltype1))
    ct2 <- paste0("density_", janitor::make_clean_names(celltype2))

    if(! ct1 %in% colnames(dens_dat))
        stop(paste0(ct1, " is not found in the data, or its density has not been computed."))
    if(! ct2 %in% colnames(dens_dat))
        stop(paste0(ct2, " is not found in the data, or its density has not been computed."))

    is.ROI <- "roi" %in% names(spe@metadata)
    if (!is.ROI) {
        dens_dat$density_coi_average <- rowMeans(as.matrix(dens_dat[, which(colnames(dens_dat) %in% c(ct1, ct2)), drop = FALSE]))
        kp <- dens_dat$density_coi_average >= quantile(dens_dat$density_coi_average, probs = probs)
        dens_dat_filter <- dens_dat[kp, ]
        rois <- data.frame(component=gl(1,sum(kp)), 
                           members=dens_dat_filter$node, 
                           x=dens_dat_filter$node_x, 
                           y=dens_dat_filter$node_y)
    } else {
        rois <- as.data.frame(spe@metadata$roi)
    }

    plotdf <- merge(rois, dens_dat,
        by.x = "members", by.y = "node",
        all.x = TRUE, sort = FALSE
    )
    plotdf <- plotdf[, c("component", c(ct1, ct2))]

    x <- rlang::sym(ct1)
    y <- rlang::sym(ct2)

    # extract aes
    aesmap <- rlang::enquos(...)
    # compute plot
    aesmap <- aesmap[!names(aesmap) %in% c("x", "y")]
    # remove x,y mappings if present

    # split aes params into those that are not aes
    # i.e. static parametrisation
    if (length(aesmap) > 0) {
        is_aes <- vapply(aesmap, rlang::quo_is_symbolic,
            FUN.VALUE = logical(1)
        )
        defaultmap <- lapply(aesmap[!is_aes], rlang::eval_tidy)
        aesmap <- aesmap[is_aes]
    } else {
        defaultmap <- list()
    }

    if (length(fit) > 1) {
        fit <- "spline"
    }

    # set some default plotting parameters
    if (is.null(defaultmap$shape)) {
        defaultmap$shape <- 21
    }

    if (is.null(defaultmap$alpha)) {
        defaultmap$alpha <- 0.8
    }

    if (isTRUE(by.roi)) {
        if (is.null(defaultmap$fill)) {
            defaultmap$fill <- "royalblue"
        }

        p <- ggplot2::ggplot(plotdf, ggplot2::aes(!!x, !!y, !!!aesmap)) +
            do.call(ggplot2::geom_point, defaultmap) +
            ggplot2::facet_wrap(~component,
                scales = "free",
                labeller = ggplot2::labeller(component = function(label) {
                    paste0("ROI #", label)
                })
            ) +
            theme_classic()
    } else {
        defaultmap$shape <- 16
        defaultmap$alpha <- 0.8

        col.p <- selectColor(length(unique(plotdf$component)))

        p <- ggplot2::ggplot(plotdf, ggplot2::aes(!!x, !!y,
            color = component, !!!aesmap
        )) +
            do.call(ggplot2::geom_point, defaultmap) +
            theme_classic() +
            scale_color_manual(name = "ROI", values = col.p)
    }

    if (fit == "spline") {
        p <- p +
            geom_smooth(
                method = "lm", formula = y ~ splines::ns(x, df = df),
                color = "red", se = FALSE
            )
    } else if (fit == "linear") {
        p <- p +
            geom_smooth(
                method = "lm", formula = y ~ x,
                color = "red", se = FALSE
            )
    } else {
        stop("The fit parameter should either be spline or linear.")
    }

    return(p)
}


utils::globalVariables(c("component"))
