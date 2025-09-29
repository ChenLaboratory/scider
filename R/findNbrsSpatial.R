#' Construct a distance-based neighbour list from cell coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param k Integer scalar for number of nearest neighbours to find. Can be used
#' with radius. See details.
#' @param radius Numeric for maximum distance to search for neighbours. Can be 
#' with k. See details
#' @param dist_func Options for distance-based weight. "idw" for inverse 
#' distance, "exp" for exponential decay, "binary" for constant weight, and 
#' "none" for raw euclidean distance.
#' @param standardisation Options for weight standardisation. "none" for 
#' nothing, and "row" for dividing weights by number of neighbours.
#' @param scale Numeric scaler for weight scaling.
#' @param nbrs_name Name of the neighbour list to be stored. Default to be "spatial".
#' @param cpu_threads Number of cpu threads for parallel computation.
#' @return A SpatialExperiment object with neighbour list stored in 
#' spe@metadata$nbrs$cell[[nbrs_name]]
#' @details
#' if only [k] is provided, neighbours are found using 
#' \link[BiocNeighbors]{findKNN}. If only [radius] is provided, neighbours are 
#' found using \link[BiocNeighbors]{findNeighbors}. If both are provided, then 
#' knn is done first then neighbours are filtered to only those within radius.
#' 
#' @export
#' @examples
#' 
#' data("xenium_bc_spe")
#' spe <- findNbrsSpatial(spe,k=20,radius=100)
findNbrsSpatial <- function(spe,
                           k = NULL,
                           radius = NULL,
                           dist_func = c("idw", "exp","binary","none"),
                           standardisation = c("none","row"),
                           scale = 1,
                           nbrs_name = NULL,
                           cpu_threads = 6) {
  if (is.null(k) && is.null(radius)) stop("k or radius is needed")
  dist_func <- match.arg(dist_func)
  standardisation <- match.arg(standardisation)
  
  get_dist <- dist_func != "binary"
  coords <- SpatialExperiment::spatialCoords(spe)
  coords <- as.matrix(coords)
  
  # Make nbrs list
  if (!is.null(k)) {# knn + radius
    knn <- BiocNeighbors::findKNN(X=coords,
                                  k=k,
                                  get.distance = get_dist || !is.null(radius),
                                  num.threads = cpu_threads)
    # Filter radius if provided. 
    n_row = nrow(coords)
    keep = `if`(is.null(knn$distance) || is.null(radius),
                matrix(TRUE,n_row,k),
                (knn$distance) <= (radius%||%Inf))
    knn$index <- knn$index[keep]
    knn$distance <- knn$distance[keep]
    
    # Convert to nbrs list
    f <- factor(row(keep)[keep],seq_len(n_row))
    nbrs <- list(index = split(knn$index,f),
                 weight = if(get_dist) split(knn$distance,f))
    # Not really needed but for consistency
    names(nbrs$index) <- NULL
    names(nbrs$weight) <- NULL
    
  } else {# radius only
    nbrs <- BiocNeighbors::findNeighbors(X=coords,
                                        threshold = radius,
                                        get.distance = get_dist,
                                        num.threads = cpu_threads)
    if (get_dist) names(nbrs)[2] <- "weight"
  }
  # Transform weight based on dist_func.
  nbrs$weight <- switch(dist_func,
                        "idw" = lapply(nbrs$weight, function(i) scale/i),
                        "exp" = lapply(nbrs$weight, function(i) exp(-i/scale)),
                        "binary" = lapply(lengths(nbrs$index), function(i) rep.int(1,i)),
                        "raw" = nbrs$weight) 
  # standardization
  nbrs$weight <- switch(standardisation,
                        "row" = lapply(nbrs$weight, function(i) i/length(i)),
                        "none" = nbrs$weight)
  
  
  if (is.null(nbrs_name)) {
    nbrs_name <- "spatial"
  }
  spe@metadata$nbrs$cell[[nbrs_name]] <- nbrs
  return(spe)
}
