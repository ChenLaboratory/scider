#' Construct a neighbour list from grid coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param n Integer. Search for neighbours within (...). Either the number of 
#' neighbors or radius
#' @param radius Numeric. Search for neighbours within the radius.
#' @param diagonal Whether to consider diagonal connection if using square grid
#' @param dist_func Options for distance-based weight. "idw" for inverse 
#' distance, "exp" for exponential decay, "binary" for constant weight, and 
#' "raw" for raw distance.
#' @param scale Numeric scaler for weight scaling.
#' @param dist_type Options of using euclidean or manhattan for distance calculation
#' @param standardisation Options for weight standardisation. "none" for 
#' nothing, and "row" for dividing weights by number of neighbours.
#' @param nbrs_name Name of the neighbour list to be stored. Default to be "grid".
#' @param cpu_threads Number of cpu threads for parallel computation.
#' @return A SpatialExperiment object with neighbour list stored in 
#' \code{spe@metadata$nbrs$grid[[nbrs_name]]}
#' 
#' @details
#' If n is used, distance is scaled to unit distance
#' @export
#' @examples
#' 
#' data("xenium_bc_spe")
#' spe <- gridDensity(spe)
#' spe <- findNbrsGrid(spe,n=3)
findNbrsGrid <- function(spe,
                         n = 1,
                         radius = NULL,
                         diagonal = FALSE,
                         dist_func = c("idw", "exp","binary","none"),
                         dist_type = c("euclidean","manhattan"),
                         standardisation = c("row","none"),
                         scale = 1,
                         nbrs_name = NULL,
                         cpu_threads = 6) {
  dist_func <- match.arg(dist_func)
  standardisation <- match.arg(standardisation)
  get_dist <- dist_func != "binary"
  if (!is.null(radius)) { # radius
    coords <- as.matrix(spe@metadata$grid_density[c("x_grid","y_grid")])
    # Fudge radius a bit due to floating point error with grid coordinates
    radius = radius*1.01
    nbrs <- BiocNeighbors::findNeighbors(X=coords,
                                         threshold = radius,
                                         get.distance = get_dist,
                                         num.threads = cpu_threads)
    if (get_dist) names(nbrs)[2] <- "weight"
  } else { # n
    coords <- as.matrix(spe@metadata$grid_density[c("node_x","node_y")])
    
    if (spe@metadata$grid_info$grid_type=="square") {
      if (diagonal) {
        nbrs_pattern = cbind(rep(seq.int(0,n),each=n),
                             rep.int(seq.int(1,n),n+1))
      } else {
        nbrs_pattern = cbind(rep.int(seq.int(0,n-1),times=seq.int(n,1)),
                             sequence(seq.int(n,1)))
      }
      nbrs_pattern = rbind(nbrs_pattern,
                           matrix(c(-nbrs_pattern[,2],nbrs_pattern[,1]),ncol=2))
      nbrs_pattern = rbind(nbrs_pattern,-nbrs_pattern)
      
      if(get_dist) {
        dist <- switch(match.arg(dist_type),
                       euclidean = sqrt(rowSums(nbrs_pattern**2)),
                       manhattan = `if`(diagonal,
                                        pmax(abs(nbrs_pattern[,1]),
                                             abs(nbrs_pattern[,2])),
                                        rowSums(abs(nbrs_pattern))))
      }
    } else {# hex
      # Convert coords to axial 
      coords[,1] = coords[,1] - (coords[,2] + bitwAnd(coords[,2],1))/2
      # Generate 1/3 of the nbrs then rotate 120 twices.
      nbrs_pattern <- matrix(c(rep(1:n,each=n+1),
                               rep.int((0:-n),times=n)),
                             ncol=2)
      nbrs_pattern <- rbind(
        nbrs_pattern,
        matrix(c(0-rowSums(nbrs_pattern),nbrs_pattern[,1]),ncol=2),
        matrix(c(nbrs_pattern[,2],0-rowSums(nbrs_pattern)),ncol=2)
      )
      
      if(get_dist) {
        dist <- switch(match.arg(dist_type),
                       euclidean = sqrt((rowSums(nbrs_pattern^2)+(0-rowSums(nbrs_pattern))^2)/2),
                       manhattan = (rowSums(abs(nbrs_pattern))+abs(rowSums(nbrs_pattern)))/2)
      }
    }
    nc = nrow(coords)
    np = nrow(nbrs_pattern)
    nbrs_all <- coords[rep(seq_len(nc),each=np),] + 
      nbrs_pattern[rep.int(seq_len(np),nc),]
    
    # Find which nbrs_all exist and which indices they should be.
    index = match(coord_hash(nbrs_all[,1],nbrs_all[,2]),
                  coord_hash(coords[,1],coords[,2]))
    keep = !is.na(index)
    f = factor(rep(seq_len(nc),each=np)[keep],seq_len(nc))
    index = index[keep]
    
    # dist 
    nbrs <- list(index = split(index,f))
    if (get_dist) {
      dist = rep.int(dist,nc)
      dist = dist[keep]
      nbrs$weight <- split(dist,f)
    }
  }
  # Transform weight based on dist_func.
  nbrs$weight <- switch(dist_func,
                        "idw" = lapply(nbrs$weight, function(i) scale/i),
                        "exp" = lapply(nbrs$weight, function(i) exp(-i/scale)),
                        "binary" = lapply(lengths(nbrs$index), function(i) rep.int(1,i)),
                        "none" = nbrs$weight)
  # standardization
  nbrs$weight <- switch(standardisation,
                        "row" = lapply(nbrs$weight, function(i) i/sum(i)),
                        "none" = nbrs$weight)
  
  # Not really needed but for consistency
  names(nbrs$index) <- NULL
  names(nbrs$weight) <- NULL
  
  nbrs_name <- nbrs_name %||% "grid"
  spe@metadata$nbrs$grid[[nbrs_name]] <- nbrs
  return(spe)
}

#' Hash two 15-bytes signed integers into one 32-bytes integer.
#' @param a,b Integer vectors of same lengths. Can be negative.
#' @details
#' Should work for a,b in range (-2^14, 2^14-1) which is good enough for our 
#' purpose. 
#' 
#' Integer in R is 32-bytes but reserve 1 byte for NA
coord_hash <- function(a,b) {
  return(bitwOr(bitwShiftL(a+16384,16),b+16384))
}