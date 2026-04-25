#' Get top highly variable genes.
#'
#' @param spe A SpatialExperiment object.
#' @param n Integer. The number of HVGs.
#' @param min.total.count Numeric. Genes with total counts less than \code{min.total.count}
#' across all cells are not considered for HVGs.
#' @param min.prop Numeric. Genes that have non-zero counts in less than the specified 
#' proportion (\code{min.prop}) of all cells are excluded from HVG selection.
#' @return A SpatialExperiment object with HVG information stored in 
#' \code{rowData(spe)$hvg} as a logical vector.
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

    # Order by mean, fit mean–dispersion trend
    ord <- order(mu)
    fit <- stats::lowess(log(mu[ord]), disp[ord])

    # Enrichment over trend
    rat <- disp[ord] / fit$y

    # Top n by dispersion ratio
    top_idx_within_sel <- head(order(rat, decreasing = TRUE), n)

    sel_idx <- which(sel)
    hvg[sel_idx[ord[top_idx_within_sel]]] <- TRUE
  }

  SummarizedExperiment::rowData(spe)$hvg <- hvg
  spe
}

