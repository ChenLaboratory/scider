#' Construct a SNN neighbour list from assay.
#'
#' @param spe A SpatialExperiment object.
#' @param assay Name of assay for clustering. Incompatible with dimred.
#' @param dimred Name of the dimensionality reduction (e.g. PCA) for clustering.
#' Incompatible with assay
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if dimred is specified.
#' @param k Integer scalar for number of nearest neighbors to find.
#' @param BNPARAM \link[BiocNeighbors]{BiocNeighborParam} object specifying the nearest neighbor algorithm. Default is Annoy.
#' @param type Type of weighting scheme for shared neighbors. Options are rank, number, and jaccard.
#' type="rank" is defined in Xu and Su (2015).
#' @param nbrs_name Name of the neighbour list to be stored in spe. Default to be 
#' assay/dimred + "_snn".
#' @param cpu_threads Number of cpu threads for parallel computation.
#' @return A spe with the clusters stored in \link[SingleCellExperiment]{reducedDims}.
#' @details
#' Construct a SNN neighbour list using either the spe's assay or reduced 
#' dimension and store it in spe@metadata$nbrs$cell
#' 
#' neighbour list contain
#' @return A SpatialExperiment object
#' @export
#' @examples
#' 
#' data("xenium_bc_spe")
#' spe <- runPCA(spe)
#' spe <- findNbrsSNN(spe,dimred="PCA")
findNbrsSNN <- function(spe,
                        assay = "counts",
                        dimred = NULL,
                        n_dimred = NULL,
                        k = 20,
                        BNPARAM = BiocNeighbors::AnnoyParam(),
                        type = c("rank", "number", "jaccard"),
                        nbrs_name = NULL,
                        cpu_threads = 6) {
  # Prep mat matrix
  if (!is.null(dimred)) {
    mat <- SingleCellExperiment::reducedDim(spe,dimred)
    if (!is.null(n_dimred)) {
      if(length(n_dimred)==1L) {
        n_dimred <- seq_len(n_dimred)
      }
      mat <- as.matrix(mat[,n_dimred,drop=FALSE])
    }
  } else {
    mat <- t(as.matrix(SummarizedExperiment::assay(spe,assay)))
  }
  
  # KNN
  print("Getting K-nearest neighbour")
  knn.args <- list(X=mat,
                     k=k,
                     num.threads=cpu_threads,
                     get.distance=FALSE,
                     BNPARAM=BNPARAM
  )
  knn <- do.call(BiocNeighbors::findKNN,knn.args)
  
  # SNN
  print("Getting shared nearest neighbour")
  type <- match.arg(type)
  snn <- .Call("C_findSNN",t(knn$index-1),k,nrow(mat),type,cpu_threads)
  snn$index <- lapply(snn$index,"+",1L)

  if (is.null(nbrs_name)) {
    if(is.null(dimred)) nbrs_name <- paste0(assay,"_snn")
    else nbrs_name <- paste0(dimred,"_snn")
  }
  spe@metadata$nbrs$cell[[nbrs_name]] <- snn
  return(spe)
}


