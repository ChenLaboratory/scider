#' Fast PCA using irlba.
#'
#' @param spe A SpatialExperiment object.
#' @param n_pcs Number of principal components to calculate
#' @param assay Name of assay used for PCA. See details for defaults.
#' @param centre Logical. Whether to centre the assay before PCA. 
#' @param scale Logical. Whether to scale the variance to 1 before PCA. 
#' @param name Name to store the PCA in the spe's \link[SingleCellExperiment]{reducedDims}
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
                   name="PCA",
                   ...) {
  if (missing(assay) && is.null(spe@assays@data[[assay]])) {
    assay <- "counts"
    message("Default assay logcounts not found. Switching to counts assay instead.")
  }
  mat <- spe@assays@data[[assay]]
  mat <- as.matrix(mat)
  n_cells = ncol(mat)
  sds = sqrt(rowSums((mat - rowMeans(mat))^2)/(n_cells-1))
  
  if (scale) {
    keep = sds!=0
    if (!all(keep)) {
      message(paste(c("Genes with 0 variance are excluded:",
                      rownames(mat)[!keep]),collapse=" "))
      }
    sds = sds[keep]
    mat = mat[keep,]
    mat = mat/sds
  }
  mat = t(mat)
  out <- irlba::irlba(mat, nv=n_pcs,
                      center=centre,
                      ...)
  pcs <- sweep(out$u, 2, out$d, "*")
  
  row.names(pcs) <- row.names(mat)
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  attr(pcs, "varExplained") <- out$d^2/(n_cells - 1)
  attr(pcs, "percentVar") <- attr(pcs, "varExplained")/sum(sds**2)*100
  SingleCellExperiment::reducedDim(spe,name) <- pcs
  return(spe)
}
