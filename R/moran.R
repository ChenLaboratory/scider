#' Calculate local Moran for 1 to 2 variables.
#'
#' @param spe A SpatialExperiment object.
#' @param data1 Numeric vector 1. Must be same length as number of grid points
#' or cells, depending on 'at'.
#' @param data2 Numeric vector 2 for bivariate local Moran. Must be same length as data1.
#' @param at Option of grid or cell for where to look for neighbour list
#' @param nbrs_name Name of the neighbour list in \code{spe@metadata$grid[[at]]}
#' for Moran's I. If NULL, will use the newest neighbour list.
#' @param hhonly Only high-high clusters, which is more interpretable. Other 
#' clusters (e.g. "Low-Low","High-Low",...) will be assigned as undefined.
#' @param significance_cutoff Cutoff for p-value to filter non-significant clusters
#' @param permutations Number of permutations for p-value.
#' @param seed Integer. For random permutations.
#' @param cpu_threads The number of cpu threads used for parallel computation.
#' 
#' @return List of lisa_value, clusters, and pseudo p-value.
#' @examples
#' data("xenium_bc_spe")
#'
#' ## At grid.
#' spe <- gridDensity(spe, coi = "Breast cancer")
#' dat <- spe@metadata$grid_density$density_breast_cancer
#' spe <- findNbrsGrid(spe)
#' res <- localMoran(spe,data1=dat, at = "grid")
#' 
#' ## At cell.
#' dat <- as.numeric(spe$cell_type=="Breast cancer")
#' spe <- findNbrsSpatial(spe,k=20)
#' res <- localMoran(spe,data1=dat,at="cell")
#' @export
#' 
localMoran <- function(spe, data1, data2=data1[],
                       at=c("grid","cell"),
                       nbrs_name = NULL,
                       hhonly = FALSE,
                       significance_cutoff = 0.05,
                       permutations = 999,
                       seed = 123456789,
                       cpu_threads=6) {
  stopifnot(length(data1)==length(data2),
            is.numeric(data1),
            is.numeric(data2))
  at <- match.arg(at)
  nbrs <- spe@metadata$nbrs[[at]][[nbrs_name %||% length(spe@metadata$nbrs[[at]])]]
  
  if (length(data1) != length(nbrs$index)) {
    stop(paste0("length of data must be same as no. of ",at,"s"))
  }
  
  r <- .Call("C_localMoran",
             nbrs$index,
             lengths(nbrs$index),
             nbrs$weight,
             data1,
             data2,
             significance_cutoff,
             permutations,
             seed,
             cpu_threads,
             hhonly)

  # Turn cluster from int to factor 
  r$cluster <- factor(r$cluster,
                      levels = 0:6,
                      labels = c("Not significant","High-High","Low-Low","Low-High",
                                 "High-Low","Undefined","Neighborless"))
  return(r)
}

#' Calculate global Moran for 1 to 2 variables.
#'
#' @param spe A SpatialExperiment object.
#' @param data1 Numeric vector 1. Must be same length as 
#' nrow(spe@metadata$grid_density) or spatialCoords(spe), depending on 'at'.
#' @param data2 Numeric vector 2 for bivariate local Moran. Must be same length as data1.
#' @param at Option of grid or cell for where to look for neighbour list
#' @param nbrs_name Name of the neighbour list in \code{spe@metadata$grid[[at]]}
#' for Moran's I
#' @param permutations Number of permutations for p-value.
#' @param seed Integer. For random permutations.
#' @param cpu_threads (optional) The number of cpu threads used for parallel
#' LISA computation
#' @return List with global lisa, p.value, and vector of permuted lisa.
#' @export
#' 
#' @examples
#' data("xenium_bc_spe")
#'
#' ## At grid.
#' spe <- gridDensity(spe, coi = "Breast cancer")
#' dat <- spe@metadata$grid_density$density_breast_cancer
#' spe <- findNbrsGrid(spe)
#' res <- globalMoran(spe,data1 = dat,at="grid")
#' res$lisa
#' ## At cell.
#' dat <- as.numeric(spe$cell_type=="Breast cancer")
#' spe <- findNbrsSpatial(spe,k=10)
#' res <- globalMoran(spe,data1 = dat,at="cell")
#' res$lisa
#' @export
globalMoran <- function(spe, data1, data2=data1[],
                        at=c("grid","cell"),
                        nbrs_name = NULL,
                        permutations = 999,
                        seed = 123456789,
                        cpu_threads=6) {
  stopifnot(length(data1)==length(data2),
            is.numeric(data1),
            is.numeric(data2))
  at <- match.arg(at)
  nbrs <- spe@metadata$nbrs[[at]][[nbrs_name %||% length(spe@metadata$nbrs[[at]])]]
  
  if (length(data1) != length(nbrs$index)) {
    stop(paste0("length of data must be same as no. of ",at,"s"))
  }
  
  r <- .Call("C_globalMoran",
             nbrs$index,
             lengths(nbrs$index),
             nbrs$weight,
             data1,
             data2,
             permutations,
             seed,
             cpu_threads)
}