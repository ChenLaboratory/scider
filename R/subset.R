#' Subset for grid level analysis
#'
#' Overwrite the default SpatialExperiment subsetting method to ensure 
#' 'grid_density' is also subsetted if 'gridLevelAnalysis' is TRUE (1 polygon 1 spot)
#' @param x A SpatialExperiment object.
#' @param i row indices for subsetting.
#' @param j col indices for subsetting.
#' @param ... further arguments to be passed to or from other methods.
#' @param drop passed on to [ indexing operator.
#' 
#' @return A SpatialExperiment object.
#'
#' @export
setMethod("[",
          c("SpatialExperiment", "ANY", "ANY"),
          function(x, i, j, ..., drop=FALSE) {
            if (missing(i)) i <- TRUE
            if (missing(j)) j <- TRUE
            x <- methods::callNextMethod()
            if ( !S4Vectors::isEmpty(SpatialExperiment::imgData(x)) ) {
              keep <- SpatialExperiment::imgData(x)$sample_id %in% unique(x$sample_id)
              SpatialExperiment::imgData(x) <- SpatialExperiment::imgData(x)[keep, ]
            }
            
            if (!is.null(x@metadata$grid_info$gridLevelAnalysis)) {
              x@metadata$grid_density = x@metadata$grid_density[j,]
            }
            return(x)
          }
)