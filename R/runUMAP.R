#' UMAP using uwot. Parameters are set to be similar to Seurat's
#'
#' @param spe A SpatialExperiment object.
#' @param n_neighbors,n_components,metric,min_dist,spread See \link[uwot]{umap}
#' @param assay Name of assay for UMAP. Incompatible with dimred.
#' @param dimred Name of the dimensionality reduction (e.g. PCA) for UMAP. Incompatible with assay
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if dimred is specified.
#' @param name Name to store the UMAP in the spe's \link[SingleCellExperiment]{reducedDims}.
#' @param ... Other parameters to be passed to \link[uwot]{umap}.
#' @return A SpatialExperiment with the UMAP stored in \link[SingleCellExperiment]{reducedDims}.
#' @export
#' @details
#' By default, runUMAP uses PCA (from \link[scider]{runPCA}). If that's 
#' unavailable, it falls back to logcounts, then counts assay. 
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- runPCA(spe)
#' spe <- runUMAP(spe,dimred="PCA",n_dimred=10)
runUMAP <- function(spe,
                    n_neighbors=30,
                    n_components=2,
                    metric="cosine",
                    min_dist=0.3,
                    spread=1,
                    assay=NULL,
                    dimred="PCA",
                    n_dimred=NULL,
                    name="UMAP",
                    ...) {
    # Getting default data for UMAP  
    if (missing(assay) && missing(dimred) && 
        (!dimred %in% SingleCellExperiment::reducedDimNames(spe))) {
        if (is.null(spe@assays@data[["logcounts"]])) {
          assay <- "counts"
        } else {
          assay <- "logcounts"
        }
        message(paste("PCA not found. Switching to",assay,"assay instead."))
    }
    
    # If both assay and dimred are provided then use dimred.
    if (!missing(assay) && missing(dimred)) {
        mat <- spe@assays@data[[assay]]
        mat <- Matrix::t(mat)
    } else {
        mat <- SingleCellExperiment::reducedDim(spe,dimred)
        if (!is.null(n_dimred)) {
            if(length(n_dimred)==1L) {
                n_dimred <- seq_len(n_dimred)
            }
            mat <- mat[,n_dimred,drop=FALSE]
        }
    }
    mat <- as.matrix(mat)
    out <- uwot::umap(X=mat,
                      n_neighbors = n_neighbors,
                      n_components = n_components,
                      metric=metric,
                      min_dist=min_dist,
                      spread=spread,
                      ...)
    colnames(out) <- paste0("UMAP", seq_len(ncol(out)))
    
    SingleCellExperiment::reducedDim(spe,name) <- out
    return(spe)
}