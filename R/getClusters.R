#' Cluster cells in spe using graph methods.
#'
#' @param spe A SpatialExperiment object.
#' @param nbrs_name Name of neighbour list for clustering. If NULL, uses the
#' newest one in spe@metadata$nbrs$cell. \link[scider]{findNbrsSNN} must have been
#' run first; otherwise getClusters errors.
#' @param method Clustering methods. Options are leiden and louvain.
#' @param resolution Higher resolution for more clusters and lower for fewer 
#' clusters. See \link[igraph]{cluster_leiden} and \link[igraph]{cluster_louvain} 
#' @param cluster_name Name to store the clusters in spe's
#' \link[SummarizedExperiment]{colData}
#' @param merge_unassigned Logical. Clusters with at most \code{min_size} cells
#' (singletons and tiny fragments, often cells left isolated by SNN pruning) are
#' set aside. If FALSE (default), their cells are labelled "unassigned". If TRUE,
#' each is merged into the retained cluster it connects to most strongly in the
#' graph, and a message reports how many cells were merged.
#' @param min_size Maximum size for a cluster to be treated as too small (see
#' merge_unassigned). Defaults to NULL, which uses either 5 or 0.01\% of the cells,
#' whichever is smaller.
#' @param seed seed for clustering
#' @param ... Other clustering arguments for \link[igraph]{cluster_leiden} or 
#' \link[igraph]{cluster_louvain} 
#' @return A spe with the clusters stored in \link[SingleCellExperiment]{reducedDims}.
#' @details
#' Cluster cells with igraph using an SNN neighbour list built by
#' \link[scider]{findNbrsSNN}, which must be run before getClusters. Any
#' neighbour list in spe@metadata$nbrs$cell can be selected via nbrs_name.
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
                        merge_unassigned = FALSE,
                        min_size = NULL,
                        seed = 1,
                        ...) {
  set.seed(seed)
  if (is.null(spe@metadata$nbrs$cell[[1]])) {
    stop("No neighbour list found. Run findNbrsSNN(spe) before getClusters().")
  }
  if (is.null(nbrs_name)){
    g <- spe@metadata$nbrs$cell[[length(spe@metadata$nbrs$cell)]]
  } else {
    g <- spe@metadata$nbrs$cell[[nbrs_name]]
  }
  g <- .nbrs2igraph(g)
  
  method <- match.arg(method)
  method.args <- list(...)
  method.args$resolution = resolution
  method.args$graph <- g
  if (method=="leiden") {
    if (is.null(method.args$objective_function))
      method.args$objective_function <- "modularity"
    # Run Leiden to convergence rather than igraph's default of 2 iterations
    if (is.null(method.args$n_iterations))
      method.args$n_iterations <- -1
  }
  
  cluster <- do.call(switch(method,
                           leiden = igraph::cluster_leiden,
                           louvain = igraph::cluster_louvain),
                    method.args)
  
  # Handle clusters too small to be meaningful
  membership <- cluster$membership
  n <- length(membership)
  count <- tabulate(membership)
  if (is.null(min_size)) min_size <- min(5, ceiling(0.0001 * n))
  is_small <- count <= min_size
  small_cells <- is_small[membership]
  n_small <- sum(small_cells)

  # Renumber the retained (big) clusters by size, high to low.
  big_ids <- which(!is_small)
  ord <- big_ids[order(count[big_ids], decreasing = TRUE)]
  relabel <- integer(length(count))
  relabel[ord] <- seq_along(ord)
  new_label <- relabel[membership]
  new_label[small_cells] <- NA_integer_

  if (n_small > 0 && merge_unassigned) {
    cells <- which(small_cells)
    new_label[cells] <- .assignNearestCluster(g, cells, new_label)
    n_merged <- n_small - sum(is.na(new_label))
    n_iso <- n_small - n_merged
    message("Merged ", n_merged, " cells from clusters with <= ", min_size,
            " cells into the nearest cluster.",
            if (n_iso > 0) paste0(" ", n_iso, " cells had no connection to a ",
                                  "retained cluster and remain 'unassigned'.") else "")
  } else if (n_small > 0) {
    message(n_small, " cells in clusters with <= ", min_size,
            " cells labelled 'unassigned'.")
  }

  labels <- ifelse(is.na(new_label), "unassigned", as.character(new_label))
  lvls <- c(as.character(seq_along(ord)), if (anyNA(new_label)) "unassigned")
  spe[[cluster_name]] <- factor(labels, levels = lvls)
  return(spe)
}


# Assign each cell in `cells` to the retained cluster it connects to most
# strongly (by summed edge weight) in graph g. `label` is a per-vertex integer
# vector of retained-cluster labels, with NA for cells not (yet) assigned.
# Cells with no edge to any retained cluster keep NA (caller treats as unassigned).
.assignNearestCluster <- function(g, cells, label) {
  inc <- igraph::incident_edges(g, cells)
  out <- rep(NA_integer_, length(cells))
  for (i in seq_along(cells)) {
    e <- inc[[i]]
    if (length(e) == 0) next
    ed <- igraph::ends(g, e, names = FALSE)
    other <- ifelse(ed[, 1] == cells[i], ed[, 2], ed[, 1])
    lab <- label[other]
    keep <- !is.na(lab)
    if (!any(keep)) next
    w <- igraph::E(g)$weight[e]
    if (is.null(w)) w <- rep(1, length(e))
    w[is.na(w)] <- 1
    agg <- tapply(w[keep], lab[keep], sum)
    out[i] <- as.integer(names(agg)[which.max(agg)])
  }
  out
}

# Convert nbrs (in spe@metadata$nbrs) into igraph's graph
# nbrs should be a list containing index & weight
.nbrs2igraph <- function(nbrs, directed=FALSE){
  interleaves <- as.vector(
    rbind(rep.int(seq_along(nbrs$index),times=lengths(nbrs$index)),
          unlist(nbrs$index)))
  # Pin the vertex count to the number of cells. Otherwise make_graph infers it
  # from max(edge id), which drops trailing cells that pruning left fully
  # isolated, making membership shorter than ncol(spe).
  g <- igraph::make_graph(interleaves, n = length(nbrs$index),
                          directed = directed) #TODO: check direction
  igraph::E(g)$weight = unlist(nbrs$weight)

  g <- igraph::simplify(g,edge.attr.comb = "first")
  return(g)
}

