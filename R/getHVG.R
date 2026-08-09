#' Get top highly variable genes.
#'
#' @param spe A SpatialExperiment object.
#' @param n Integer. The number of HVGs.
#' @param min.total.count Numeric. Genes with total counts less than \code{min.total.count}
#' across all cells are not considered for HVGs.
#' @param min.prop Numeric. Genes that have non-zero counts in less than the specified
#' proportion (\code{min.prop}) of all cells are excluded from HVG selection.
#' @return A SpatialExperiment object with HVG information stored in
#' \code{rowData(spe)$hvg} as a logical vector. When a mean-dispersion trend is
#' fitted (i.e. more eligible genes than requested HVGs), the per-gene trend is
#' also stored in \code{metadata(spe)$hvg_trend} for use with \code{plotHVG}.
#'
#' @details
#' \code{getHVG} adopts a fast approach of NB dispersion estimation for all the genes across
#' all cells. A lowess curve is fit to represent mean-dispersion trend. Top HVGs are selected
#' based on the ratio of gene-wise dispersion and their trended dispersion.
#'
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- getHVG(spe, n=100)
getHVG <- function(spe,
                   n = 1000,
                   min.total.count = 100,
                   min.prop = 0.01) {
  ngenes <- nrow(spe)
  if (n >= ngenes) {
    message("n >= the total number of genes. All genes are used as HVGs.")
    SummarizedExperiment::rowData(spe)$hvg <- rep(TRUE, ngenes)
    return(spe)
  }

  # Initialise HVG flag
  hvg <- logical(ngenes)

  # Get count matrix once
  cnt <- spe@assays@data$counts

  # Filter out lowly expressed genes
  total_count <- Matrix::rowSums(cnt)
  prop_nz     <- Matrix::rowSums(cnt > 0) / ncol(cnt)

  sel <- total_count >= min.total.count & prop_nz >= min.prop
  n_sel <- sum(sel)

  if (n_sel == 0L) {
    message("No genes pass the specified thresholds. No HVGs selected.")
    SummarizedExperiment::rowData(spe)$hvg <- hvg
    return(spe)
  }

  if (n >= n_sel) {
    message(n_sel, " genes pass the specified thresholds. ",
            "All these genes are used as HVGs.")
    hvg[sel] <- TRUE
  } else {
    mat <- cnt[sel, , drop = FALSE]

    # mean & variance per gene
    mu  <- Matrix::rowMeans(mat)
    s2  <- Matrix::rowMeans(mat^2) - mu^2

    # NB method-of-moments dispersion: disp = (s2 - mu) / mu^2
    disp <- (s2 - mu) / pmax(mu^2, 1e-8)
    disp[disp < 1e-4] <- 1e-4  # lower bound for stability

    # Order by mean, fit mean-dispersion trend
    ord <- order(mu)
    fit <- stats::lowess(log(mu[ord]), disp[ord])

    # Enrichment over trend
    rat <- disp[ord] / fit$y

    # Top n by dispersion ratio
    top_idx_within_sel <- head(order(rat, decreasing = TRUE), n)

    sel_idx <- which(sel)
    hvg[sel_idx[ord[top_idx_within_sel]]] <- TRUE

    # Store the per-gene mean-dispersion trend (eligible genes only) so it can
    # be visualised with plotHVG().
    trend <- numeric(n_sel)
    trend[ord] <- fit$y
    hvg_sel <- logical(n_sel)
    hvg_sel[ord[top_idx_within_sel]] <- TRUE
    spe@metadata$hvg_trend <- data.frame(
      mean_expr  = as.numeric(mu),
      dispersion = as.numeric(disp),
      trend      = trend,
      hvg        = hvg_sel,
      row.names  = rownames(mat)
    )
  }

  SummarizedExperiment::rowData(spe)$hvg <- hvg
  spe
}


#' Plot the mean-dispersion trend from getHVG().
#'
#' Visualises each eligible gene's dispersion against its mean expression,
#' overlaid with the fitted lowess trend, and highlights the selected highly
#' variable genes. Requires \code{\link{getHVG}} to have been run with a fitted
#' trend (i.e. more eligible genes than requested HVGs).
#'
#' @param spe A SpatialExperiment processed by \code{\link{getHVG}}.
#' @param pt.size Point size. Default 0.6.
#' @param pt.alpha Point alpha (0-1). Default 0.6.
#' @param cols Length-2 vector of colours for non-HVG and HVG points. Default
#' \code{c("grey70", "firebrick")}.
#' @param line.col Colour of the fitted trend line. Default "blue".
#' @return A ggplot object (dispersion vs mean expression, both on log10 axes).
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- getHVG(spe, n = 100)
#' plotHVG(spe)
plotHVG <- function(spe,
                    pt.size = 0.6,
                    pt.alpha = 0.6,
                    cols = c("grey70", "firebrick"),
                    line.col = "blue") {
  df <- spe@metadata$hvg_trend
  if (is.null(df)) {
    stop("No mean-dispersion trend found in metadata(spe)$hvg_trend. Run ",
         "getHVG() with a fitted trend (n smaller than the number of eligible ",
         "genes) first.")
  }
  df <- df[order(df$mean_expr), ]
  ggplot2::ggplot(df, ggplot2::aes(x = mean_expr, y = dispersion)) +
    ggplot2::geom_point(ggplot2::aes(color = hvg), size = pt.size,
                        alpha = pt.alpha) +
    ggplot2::geom_line(ggplot2::aes(y = trend), color = line.col,
                       linewidth = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(
      values = stats::setNames(cols, c("FALSE", "TRUE")),
      labels = c("Other", "HVG"), name = NULL) +
    ggplot2::labs(x = "Mean expression", y = "Dispersion") +
    ggplot2::theme_classic()
}

utils::globalVariables(c("mean_expr", "dispersion", "trend", "hvg"))

