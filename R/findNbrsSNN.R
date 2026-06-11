#' Construct a SNN neighbour list from assay.
#'
#' @param spe A SpatialExperiment object.
#' @param assay Name of assay for clustering. Incompatible with dimred.
#' @param dimred Name of the dimensionality reduction (e.g. PCA) for clustering.
#' Incompatible with assay
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' dimred is specified. Defaults to NULL, which uses all available dimensions of
#' the reduced dim (e.g. all PCs produced by runPCA).
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
                        assay = NULL,
                        dimred = "PCA",
                        n_dimred = NULL,
                        k = 20,
                        BNPARAM = BiocNeighbors::AnnoyParam(),
                        type = c("rank", "number", "jaccard"),
                        nbrs_name = NULL,
                        cpu_threads = 6) {
  # Getting default data  
  if (missing(assay) && missing(dimred) && 
      (!dimred %in% SingleCellExperiment::reducedDimNames(spe))) {
    dimred=NULL
    if (is.null(spe@assays@data[["logcounts"]])) {
      assay <- "counts"
    } else {
      assay <- "logcounts"
    }
    warning(paste("PCA not found. Switching to",assay,"assay instead."))
  }

  if (!is.null(dimred)) {
    mat <- SingleCellExperiment::reducedDim(spe,dimred)
    if (is.null(n_dimred)) n_dimred <- ncol(mat)
    if (length(n_dimred)==1L) {
      if (n_dimred > ncol(mat)) {
        message("n_dimred (", n_dimred, ") exceeds available dimensions (",
                ncol(mat), "). Using all ", ncol(mat), ".")
        n_dimred <- ncol(mat)
      }
      n_dimred <- seq_len(n_dimred)
    }
    mat <- mat[,n_dimred,drop=FALSE]
  } else {
    mat <- Matrix::t(SummarizedExperiment::assay(spe,assay))
  }
  
  # KNN
  if (k > ncol(spe)-1) {
    message(paste("k is larger than number of points. Setting k to",ncol(spe)-1))
    k = ncol(spe)-1
  }
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


