#' Fast PCA using irlba.
#'
#' @param spe A SpatialExperiment object.
#' @param n_pcs Number of principal components to calculate
#' @param assay Name of assay used for PCA. See details for defaults.
#' @param centre Logical. Whether to centre the assay before PCA.
#' @param scale Logical. Whether to scale the variance to 1 before PCA.
#' @param clip Maximum absolute z-score after scaling. Values beyond this are
#' clipped. Prevents rare-marker genes (near-zero SD) from creating extreme
#' outliers that fragment the UMAP. Defaults to \code{sqrt(ncol(spe))}. Set to 
#' \code{Inf} to disable.
#' @param name Name to store the PCA in the spe's \link[SingleCellExperiment]{reducedDims}
#' @param genes Subset of features for PCA. Can be a column in rowData or a vector 
#' of gene names, indices, or booleans. Default to hvg if \link[scider]{getHVG} was run.
#' @param ... Other parameters to be passed to \link[irlba]{irlba}.
#' @return A SpatialExperiment with the PCA stored in \link[SingleCellExperiment]{reducedDims}.
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe <- runPCA(spe)
#' @details
#' By default, runPCA uses logcounts assay (from \link[scider]{normalizeAssay}). 
#' If that's unavailable, it falls back to counts assay
#' 
#' @export
runPCA <- function(spe,
                   n_pcs=50,
                   assay="logcounts",
                   centre = TRUE,
                   scale = TRUE,
                   clip = NULL,
                   name="PCA",
                   genes="hvg",
                   ...) {
  if (missing(assay) && is.null(spe@assays@data[[assay]])) {
    assay <- "counts"
    message("Default assay logcounts not found. Switching to counts assay instead.")
  }
  mat <- spe@assays@data[[assay]]
  
  # Subset to only the relevant genes.
  if (length(genes)==1 && is.character(genes)) {
    genes <- SummarizedExperiment::rowData(spe)[[genes]]
  }
  if (!is.null(genes)) mat <- mat[genes,,drop=FALSE]
    
  n_cells <- ncol(mat)

  mu  <- Matrix::rowMeans(mat)
  ex2 <- Matrix::rowMeans(mat^2)

  var_hat <- (ex2 - mu^2) * n_cells / (n_cells - 1)
  var_hat[var_hat < 0] <- 0  # numerical guard
  sds <- sqrt(var_hat)
  
  if (scale) {
    keep <- sds != 0 & !is.na(sds)
    if (!all(keep)) {
      message(paste(c("Genes with 0 variance are excluded:",
                      rownames(mat)[!keep]),collapse=" "))
    }
    mat <- mat[keep, , drop=FALSE]
    sds <- sds[keep]

    # Row-wise scaling; works for dense and sparse matrices
    mat <- mat / sds

    # Clip extreme z-scores.
    if (is.null(clip)) clip <- sqrt(n_cells)
    mat@x <- pmin(pmax(mat@x, -clip), clip)

    # After scaling, all kept genes have unit variance
    sds <- rep.int(1, nrow(mat))
  }

  # irlba expects samples in rows -> transpose
  mat <- Matrix::t(mat)

  n_pcs <- min(n_pcs, min(dim(mat)) - 1L)
  if (n_pcs > floor(ncol(mat) * 0.3))
    message("n_pcs (", n_pcs, ") exceeds 30% of the number of genes (", ncol(mat),
            "). irlba approximation quality may degrade for higher components.",
            " Consider reducing n_pcs.")

  out <- irlba::irlba(mat, nv = n_pcs, center = centre, ...)

  pcs <- sweep(out$u, 2, out$d, "*")
  row.names(pcs) <- row.names(mat)
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))

  varExplained <- out$d^2 / (n_cells - 1)
  attr(pcs, "varExplained") <- varExplained
  attr(pcs, "percentVar")   <- varExplained / sum(sds^2) * 100

  SingleCellExperiment::reducedDim(spe, name) <- pcs
  return(spe)
}
