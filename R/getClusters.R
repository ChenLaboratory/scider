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
#' @param unassigned How to handle cells in clusters with at most \code{min_size}
#' cells (singletons and tiny fragments, often left isolated by SNN pruning). One
#' of: "merge" (default) graph-merges each cell into its most-connected retained
#' cluster; "label" leaves all such cells "unassigned"; "discard" removes them 
#' from the spe (changing ncol(spe)).
#' @param min_size Maximum size for a cluster to be treated as too small (see
#' unassigned). Defaults to NULL, which uses either 5 or 0.01\% of the cells,
#' whichever is smaller.
#' @param start_from Integer at which cluster numbering starts. 1 (default)
#' numbers clusters 1..K.
#' @param seed seed for clustering
#' @param verbose Logical. Whether to report how cells in tiny clusters were
#' reassigned/labelled/discarded (see unassigned). Defaults to FALSE.
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
                        unassigned = c("merge", "label", "discard"),
                        min_size = NULL,
                        start_from = 1,
                        seed = 1,
                        verbose = FALSE,
                        ...) {
  set.seed(seed)
  unassigned <- match.arg(unassigned)
  if (is.null(spe@metadata$nbrs$cell[[1]])) {
    stop("No neighbour list found. Run findNbrsSNN(spe) before getClusters().")
  }
  if (is.null(nbrs_name)){
    nbrs <- spe@metadata$nbrs$cell[[length(spe@metadata$nbrs$cell)]]
  } else {
    nbrs <- spe@metadata$nbrs$cell[[nbrs_name]]
  }
  g <- .nbrs2igraph(nbrs)
  
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

  # Renumber the retained (big) clusters by size, high to low. Numbering starts
  # at start_from (1 by default).
  big_ids <- which(!is_small)
  ord <- big_ids[order(count[big_ids], decreasing = TRUE)]
  cluster_ids <- seq_along(ord) - 1L + start_from
  relabel <- integer(length(count))
  relabel[ord] <- cluster_ids
  new_label <- relabel[membership]
  new_label[small_cells] <- NA_integer_

  if (n_small > 0) {
    switch(unassigned,
      merge = {
        cells <- which(small_cells)
        new_label[cells] <- .assignNearestCluster(g, cells, new_label)
        n_graph <- n_small - sum(is.na(new_label))
        leftover <- which(is.na(new_label))
        if (length(leftover)) {
          coords <- .nbrsSpace(spe, nbrs)
          new_label[leftover] <-
            .assignNearestCentroid(coords, new_label, leftover)
        }
        if (verbose)
          message("Reassigned all ", n_small, " cells from clusters with <= ",
                min_size, " cells (", n_graph, " by graph, ", length(leftover),
                " by distance).")
      },
      label = {
        if (verbose)
          message(n_small, " cells in clusters with <= ", min_size,
                " cells labelled 'unassigned'.")
      },
      discard = {
        if (verbose)
          message("Discarded ", n_small, " cells in clusters with <= ", min_size,
                " cells. Neighbour list cleared; re-run findNbrsSNN() to ",
                "re-cluster the remaining cells.")
      })
  }

  if (unassigned == "discard" && n_small > 0) {
    keep_cells <- !small_cells
    spe <- spe[, keep_cells]
    new_label <- new_label[keep_cells]
    # The stored cell neighbour lists no longer match the subset cells; drop them
    # so a later getClusters errors (forcing a fresh findNbrsSNN) instead of
    # silently using a stale graph.
    spe@metadata$nbrs$cell <- NULL
  }

  labels <- ifelse(is.na(new_label), "unassigned", as.character(new_label))
  lvls <- c(as.character(cluster_ids), if (anyNA(new_label)) "unassigned")
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

# Assign each cell in `cells` to the nearest retained-cluster centroid in the
# coordinate space `coords` (cells x dims). `label` is a per-cell integer vector
# of retained-cluster labels (NA for cells not yet assigned). Centroids are the
# column means of the currently-assigned cells in each cluster.
.assignNearestCentroid <- function(coords, label, cells) {
  assigned <- !is.na(label)
  cl <- sort(unique(label[assigned]))
  cent <- vapply(cl, function(k)
    colMeans(coords[assigned & label == k, , drop = FALSE]),
    numeric(ncol(coords)))            # dims x n_clusters
  out <- integer(length(cells))
  for (i in seq_along(cells)) {
    d2 <- colSums((cent - coords[cells[i], ])^2)
    out[i] <- cl[which.min(d2)]
  }
  out
}

# Retrieve the coordinate matrix (cells x dims) that built the neighbour graph,
# using the dimred/assay recorded by findNbrsSNN.
.nbrsSpace <- function(spe, nbrs) {
  if (!is.null(nbrs$dimred)) {
    m <- as.matrix(SingleCellExperiment::reducedDim(spe, nbrs$dimred))
    if (!is.null(nbrs$dims)) m <- m[, nbrs$dims, drop = FALSE]
    m
  } else if (!is.null(nbrs$assay)) {
    as.matrix(Matrix::t(SummarizedExperiment::assay(spe, nbrs$assay)))
  } else {
    stop("The coordinate space for the neighbour graph was not recorded; ",
         "re-run findNbrsSNN() so unassigned=\"force\" can use it.")
  }
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

