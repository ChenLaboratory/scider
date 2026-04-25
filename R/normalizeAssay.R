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
                           scale.factor = 1e4,
                           assay = "counts",
                           name = "logcounts") {
  assay <- SummarizedExperiment::assay(spe,assay)
  library_size <- scale.factor/Matrix::colSums(assay)
  assay <- Matrix::colScale(assay,library_size)
  switch(match.arg(transformation),
         log = {
           assay <- log1p(assay)
         })
  SummarizedExperiment::assay(spe,name) <- assay
  return(spe)
}