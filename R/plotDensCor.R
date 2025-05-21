#' Plot density correlation between two cell types
#'
#' @param spe A SpatialExperiment object.
#' @param celltype1 Cell type 1 to compare.
#' @param celltype2 Cell type 2 to compare.
#' @param roi Character. The name of the group or cell type on which
#' the roi is computed. Default is NULL for no facetting by ROI 
#' @param probs A numeric scalar. The threshold of proportion that used to
#' filter grids by density when ROIs have not been identified previously.
#' Ignored if 'roi' is present in the 'metadata' component of spe. Default to 0.85.
#' @param fit Character. Options are "spline" and "linear".
#' @param df Integer. Degrees of freedom of the spline fit.
#' Default to 3 (i.e., a cubic spline fit).
#' @param pt.shape shape of points.
#' @param pt.size size of points.
#' @param pt.alpha alpha of points between 0 and 1.
#' @param line.type shape of line.
#' @param line.width size of line.
#' @param line.alpha alpha of line between 0 and 1.
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
#' plotDensCor(spe, celltype1 = "Breast cancer", celltype2 = "Fibroblasts", roi = coi)
#'
plotDensCor <- function(spe, celltype1 = NULL, celltype2 = NULL,
                        roi = NULL, probs = 0.85,
                        fit = c("spline", "linear"), df = 3, 
                        pt.shape = 21, 
                        pt.size = 1.5, 
                        pt.alpha = 1,
                        line.type = 1,
                        line.width = 1,
                        line.alpha = 1) {
    if (!("grid_density" %in% names(spe@metadata))) {
        stop("Please run gridDensity before using this function.")
    }
    fit <- match.arg(fit)
    
    dens_dat <- as.data.frame(spe@metadata$grid_density)
    
    # clean names
    ct1 <- paste0("density_", janitor::make_clean_names(celltype1))
    ct2 <- paste0("density_", janitor::make_clean_names(celltype2))
    
    if(! ct1 %in% colnames(dens_dat))
        stop(paste0(ct1, " is not found in the data, or its density has not been computed."))
    if(! ct2 %in% colnames(dens_dat))
        stop(paste0(ct2, " is not found in the data, or its density has not been computed."))
    
    if (is.null(roi)) {
        dens_dat$density_coi_average <- (dens_dat[[ct1]]+dens_dat[[ct2]])/2
        kp <- dens_dat$density_coi_average >= quantile(dens_dat$density_coi_average, probs = probs)
        dens_dat_filter <- dens_dat[kp, ]
        rois <- data.frame(component=gl(1,sum(kp)), 
                           members=dens_dat_filter$node)
    } else {
        roi <- gsub("_roi$", "", roi)
        roi <- janitor::make_clean_names(roi)
        roi <- paste(c(sort(roi),"roi"), collapse="_")
        if (is.null(spe@metadata[[roi]])) {
            stop(paste(
                roi, " is not found in metadata of spe. Please run
                findROI() first."
            ))
        }
        rois <- as.data.frame(spe@metadata[[roi]])
    }
    
    plotdf <- merge(rois, dens_dat,
                    by.x = "members", by.y = "node",
                    all.x = TRUE, sort = FALSE
    )
    plotdf <- plotdf[, c("component", c(ct1, ct2))]
    
    p <- ggplot2::ggplot(plotdf, ggplot2::aes(.data[[ct1]], .data[[ct2]])) +
        ggplot2::geom_point(shape = pt.shape,
                            size = pt.size,
                            alpha = pt.alpha,
                            fill = "royalblue") +
        geom_line(stat="smooth", method = "lm",
                  formula = switch(fit,
                                   spline = y ~ splines::ns(x, df = df),
                                   linear = y ~ x),
                  color = "red", se = FALSE,
                  linetype=line.type,
                  linewidth=line.width,
                  alpha=line.alpha) +
        theme_classic()
    if (!is.null(roi)) {
        p <- p + ggplot2::facet_wrap(~component,
                                    scales = "free",
                                    labeller = ggplot2::labeller(component = function(label) {
                                        paste0("ROI #", label)
                                    }))
    }
    return(p)
}
utils::globalVariables(c("component"))
