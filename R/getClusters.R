#' Cluster cells in spe using graph methods.
#'
#' @param spe A SpatialExperiment object.
#' @param nbrs_name Name of neighbour list for clustering. If NULL, will use 
#' the newest one in spe@metadata$nbrs$cell or create one if none are available.
#' @param method Clustering methods. Options are leiden and louvain.
#' @param resolution Higher resolution for more clusters and lower for fewer 
#' clusters. See \link[igraph]{cluster_leiden} and \link[igraph]{cluster_louvain} 
#' @param cluster_name Name to store the clusters in spe's 
#' \link[SummarizedExperiment]{colData}
#' @param seed seed for clustering
#' @param ... Other clustering arguments for \link[igraph]{cluster_leiden} or 
#' \link[igraph]{cluster_louvain} 
#' @return A spe with the clusters stored in \link[SingleCellExperiment]{reducedDims}.
#' @details
#' Cluster cells with igraph using SNN calculated by \link[scider]{findNbrsSNN}.
#' Any neighbour list in spe@metadata$nbrs$cell can also be used
#' @return A SpatialExperiment object
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- normalizeAssay(spe)
#' spe <- runPCA(spe)
#' spe <- findNbrsSNN(spe,dimred="PCA")
#' spe <- getClusters(spe, resolution=0.5)

getClusters <- function(spe,
                        nbrs_name = NULL,
                        method = c("leiden", "louvain"),
                        resolution = 1,
                        cluster_name = "cluster",
                        seed = 1,
                        ...) {
  set.seed(seed)
  if (is.null(nbrs_name)){
    if (is.null(spe@metadata$nbrs$cell[[1]])) {
      message("No neighbour list found. Calculating snn.")
      spe <- findNbrsSNN(spe)
    }
    g <- spe@metadata$nbrs$cell[[length(spe@metadata$nbrs$cell)]]
  } else {
    g <- spe@metadata$nbrs$cell[[nbrs_name]]
  }
  g <- .nbrs2igraph(g)
  
  method <- match.arg(method)
  method.args <- list(...)
  method.args$resolution = resolution
  method.args$graph <- g
  if (method=="leiden" && is.null(method.args$objective_function)) {
    method.args$objective_function <- "modularity"
  }
  
  cluster <- do.call(switch(method,
                           leiden = igraph::cluster_leiden,
                           louvain = igraph::cluster_louvain),
                    method.args)
  
  # Rename clusters based on their sizes high to low
  count <- tabulate(cluster$membership)
  map <- order(order(count,decreasing=TRUE))
  spe[[cluster_name]] <- factor(map[cluster$membership])
  return(spe)
}


# Convert nbrs (in spe@metadata$nbrs) into igraph's graph 
# nbrs should be a list containing index & weight
.nbrs2igraph <- function(nbrs, directed=FALSE){
  interleaves <- as.vector(
    rbind(rep.int(seq_along(nbrs$index),times=lengths(nbrs$index)),
          unlist(nbrs$index)))
  g <- igraph::make_graph(interleaves,directed=directed)
  igraph::E(g)$weight = unlist(nbrs$weight)

  g <- igraph::simplify(g,edge.attr.comb = "first")
  return(g)
}

