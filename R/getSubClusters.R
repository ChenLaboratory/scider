#' Sub-cluster a single cluster.
#'
#' Take one cluster from an existing clustering and partition it further by
#' community detection on the sub-graph of its cells, similar to
#' \code{Seurat::FindSubCluster}. New sub-clusters are labelled \code{<cluster>_1},
#' \code{<cluster>_2}, ... (largest first).
#'
#' @param spe A SpatialExperiment object.
#' @param cluster The cluster label to sub-cluster (e.g. 2 or "2").
#' @param resolution Resolution for the sub-clustering. Higher gives more
#' sub-clusters. See \link[igraph]{cluster_leiden} and \link[igraph]{cluster_louvain}.
#' @param nbrs_name Name of neighbour list to use. If NULL, uses the newest one in
#' spe@metadata$nbrs$cell. \link[scider]{findNbrsSNN} must have been run first.
#' @param method Clustering method. Options are leiden and louvain.
#' @param cluster_name Name of the cluster column in
#' \link[SummarizedExperiment]{colData} to sub-cluster.
#' @param new_name Name of the column to store the result. Defaults to
#' overwriting cluster_name.
#' @param relabel Logical. If TRUE, renumber all clusters and sub-clusters by size
#' (largest first) after sub-clustering, replacing the \code{<cluster>_<n>} labels
#' with plain numbers. "unassigned" is kept as the last level. Defaults to FALSE.
#' @param min_size Sub-clusters with at most this many cells are merged into the
#' larger sub-clusters (graph connection first, then nearest centroid), as in
#' \link[scider]{getClusters}. Defaults to NULL, which uses either 5 or 0.01\% of
#' the cells, whichever is smaller.
#' @param start_from Integer at which numbering starts, for both the sub-cluster
#' suffixes and the relabel renumbering. Defaults to NULL, which infers the base
#' (0- or 1-based) from the existing labels.
#' @param sep Separator between a cluster and its sub-cluster number. Default "_".
#' @param seed Seed for clustering.
#' @param ... Other clustering arguments for \link[igraph]{cluster_leiden} or
#' \link[igraph]{cluster_louvain}.
#' @return A SpatialExperiment with the sub-clusters stored in colData.
#' @details
#' The sub-graph of the target cluster's cells is induced from the SNN graph and
#' clustered on its own. Cells that are weakly connected within the cluster may
#' form small or singleton sub-clusters; lower the resolution if this is not
#' wanted.
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- normalizeAssay(spe)
#' spe <- runPCA(spe)
#' spe <- findNbrsSNN(spe, dimred = "PCA")
#' spe <- getClusters(spe, resolution = 0.5)
#' # Split cluster 2 into 2_1, 2_2, ...
#' spe <- getSubClusters(spe, cluster = 2, resolution = 0.5)
getSubClusters <- function(spe,
                           cluster,
                           resolution = 1,
                           nbrs_name = NULL,
                           method = c("leiden", "louvain"),
                           cluster_name = "cluster",
                           new_name = cluster_name,
                           relabel = FALSE,
                           min_size = NULL,
                           start_from = NULL,
                           sep = "_",
                           seed = 1,
                           ...) {
  set.seed(seed)
  method <- match.arg(method)
  if (!cluster_name %in% colnames(SummarizedExperiment::colData(spe)))
    stop("'", cluster_name, "' not found in colData(spe).")
  if (is.null(spe@metadata$nbrs$cell[[1]]))
    stop("No neighbour list found. Run findNbrsSNN(spe) before getSubClusters().")

  cl <- spe[[cluster_name]]
  cl_chr <- as.character(cl)
  target <- as.character(cluster)
  cells <- which(cl_chr == target)
  if (length(cells) < 2L)
    stop("Cluster '", target, "' has fewer than 2 cells to sub-cluster.")

  # Numbering base for sub-cluster suffixes and relabel; infer if not given.
  if (is.null(start_from)) {
    orig_num <- suppressWarnings(as.numeric(setdiff(unique(cl_chr), "unassigned")))
    start_from <- if (any(!is.na(orig_num))) min(orig_num, na.rm = TRUE) else 1
  }

  # Induce and cluster the sub-graph of the target cluster's cells.
  if (is.null(nbrs_name)) {
    nbrs <- spe@metadata$nbrs$cell[[length(spe@metadata$nbrs$cell)]]
  } else {
    nbrs <- spe@metadata$nbrs$cell[[nbrs_name]]
  }
  sub_g <- igraph::induced_subgraph(.nbrs2igraph(nbrs), cells)

  method.args <- list(...)
  method.args$resolution <- resolution
  method.args$graph <- sub_g
  if (method == "leiden") {
    if (is.null(method.args$objective_function))
      method.args$objective_function <- "modularity"
    if (is.null(method.args$n_iterations))
      method.args$n_iterations <- -1
  }
  sub <- do.call(switch(method, leiden = igraph::cluster_leiden,
                                louvain = igraph::cluster_louvain), method.args)

  # Number retained sub-clusters by size (largest first), suffix from start_from.
  # Small sub-clusters (<= min_size) are merged in, matching
  # getClusters(unassigned = "merge"): graph connection first, then nearest
  # centroid for any cell the subgraph could not place.
  sub_mem <- sub$membership
  sub_count <- tabulate(sub_mem)
  if (is.null(min_size)) min_size <- min(5, ceiling(0.0001 * ncol(spe)))
  is_small_sub <- sub_count <= min_size
  big_sub <- which(!is_small_sub)
  if (length(big_sub) == 0) {
    warning("Cluster '", target, "' produced only sub-clusters of <= ", min_size,
            " cells; leaving it unchanged. Try a lower resolution.")
    spe[[new_name]] <- cl
    return(spe)
  }
  sub_ord <- big_sub[order(sub_count[big_sub], decreasing = TRUE)]
  sub_relabel <- integer(length(sub_count))
  sub_relabel[sub_ord] <- seq_along(sub_ord) - 1L + start_from
  sub_id <- sub_relabel[sub_mem]
  small_sub_cells <- is_small_sub[sub_mem]
  sub_id[small_sub_cells] <- NA_integer_

  n_small <- sum(small_sub_cells)
  if (n_small > 0) {
    sc <- which(small_sub_cells)
    sub_id[sc] <- .assignNearestCluster(sub_g, sc, sub_id)
    leftover <- which(is.na(sub_id))
    if (length(leftover)) {
      coords <- .nbrsSpace(spe, nbrs)[cells, , drop = FALSE]
      sub_id[leftover] <- .assignNearestCentroid(coords, sub_id, leftover)
    }
  }

  new_chr <- cl_chr
  new_chr[cells] <- paste(target, sub_id, sep = sep)

  message("Split cluster '", target, "' into ", length(sub_ord), " sub-clusters",
          if (n_small > 0) paste0(" (", n_small, " cells from sub-clusters <= ",
                                  min_size, " merged)") else "", ".")

  if (relabel) {
    # Renumber all clusters by size, preserving the numbering base.
    keep <- new_chr != "unassigned"
    ord  <- names(sort(table(new_chr[keep]), decreasing = TRUE))
    ids  <- seq_along(ord) - 1L + start_from
    map  <- ids
    names(map) <- ord
    new_chr[keep] <- as.character(map[new_chr[keep]])
    lvls <- c(as.character(ids), if (any(!keep)) "unassigned")
  } else {
    # Keep the existing level order, replacing the target cluster with its
    # sub-cluster labels (in sub-id order).
    orig_levels <- if (is.factor(cl)) levels(cl) else .orderClusterLevels(cl_chr)
    sub_lab_sorted <- paste(target, sort(unique(sub_id)), sep = sep)
    pos <- match(target, orig_levels)
    lvls <- append(orig_levels[-pos], sub_lab_sorted, after = pos - 1L)
  }

  spe[[new_name]] <- factor(new_chr, levels = lvls)
  return(spe)
}
