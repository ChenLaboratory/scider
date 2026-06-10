#' Perform log normalization for counts
#'
#' @param spe A SpatialExperiment object.
#' @param transformation Choice of transformation. "Log" for log1p
#' @param scale.factor Factor to multiply the count of each cell by. A single 
#' value or a numeric vector of length equal to the number of cells.
#' @param assay Name of assay in spe to perform the transformation on
#' @param name Name of the transformed assay
#' @return A SpatialExperiment object 
#' @export
#' @examples
#' data("xenium_bc_spe")
#' spe <- normalizeAssay(spe)
normalizeAssay <- function(spe,
                           transformation = c("log"),
                           scale.factor = 1e2,
                           assay = "counts",
                           name = "logcounts") {
  mat <- SummarizedExperiment::assay(spe, assay)
  dn <- dimnames(mat)
  library_size <- scale.factor / Matrix::colSums(mat)
  # Multiply each column by its per-cell scaling factor.
  # Using %*% Diagonal keeps the result as a sparse dgCMatrix and avoids
  # Matrix::colScale, which was removed in Matrix >= 1.6 and returned a dense
  # dgeMatrix in earlier versions (causing a 'Dimnames' S4 slot error on
  # downstream assay assignment).
  # dimnames must be restored afterwards because %*% drops them.
  mat <- mat %*% Matrix::Diagonal(x = library_size)
  mat <- methods::as(mat, "dgCMatrix")
  dimnames(mat) <- dn
  switch(match.arg(transformation),
         log = {
           mat <- log1p(mat)
         })
  SummarizedExperiment::assay(spe, name) <- mat
  return(spe)
}
