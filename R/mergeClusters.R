#' Manually merge clusters.
#'
#' Combine specified clusters into a single cluster, optionally relabelling all
#' clusters afterwards. This is a manual alternative to the modularity-based
#' \link[bluster]{mergeCommunities}: the user chooses exactly which clusters to
#' combine.
#'
#' @param spe A SpatialExperiment object.
#' @param merge A vector of cluster labels to merge into one (e.g. c(2, 5, 7)),
#' or a list of such vectors to perform several merges in one call.
#' @param cluster_name Name of the cluster column in
#' \link[SummarizedExperiment]{colData} to merge.
#' @param new_name Name of the column to store the result. Defaults to
#' overwriting cluster_name.
#' @param merge_label Optional label for the merged cluster (single merge only).
#' Defaults to the merged cluster labels joined by \code{sep}, e.g. "2&5&7".
#' @param relabel Logical. If TRUE, renumber all clusters by size (largest first)
#' after merging, replacing the merged label with a number. "unassigned" is always
#' kept as the last level. Defaults to FALSE (keep the merged label).
#' @param start_from Integer at which renumbering starts when relabel=TRUE.
#' Defaults to NULL, which infers the base (0- or 1-based) from the existing
#' labels, matching what getClusters produced. Ignored when relabel=FALSE.
#' @param sep Separator for the auto-generated merged label. Default "&".
#' @return A SpatialExperiment with the merged clusters stored in colData.
#' @details
#' Clusters listed in \code{merge} that are not present in cluster_name are
#' ignored; a merge group with fewer than two present clusters is skipped with a
#' warning. Any "unassigned" cells are left untouched (and kept as the last
#' factor level) unless explicitly included in \code{merge}.
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- normalizeAssay(spe)
#' spe <- runPCA(spe)
#' spe <- findNbrsSNN(spe, dimred = "PCA")
#' spe <- getClusters(spe, resolution = 0.5)
#' # Merge clusters 2, 5 and 7 into one labelled "2&5&7"
#' spe <- mergeClusters(spe, merge = c(2, 5, 7))
#' # Merge and then renumber everything 1..K by size
#' spe <- mergeClusters(spe, merge = c(2, 5, 7), relabel = TRUE)
mergeClusters <- function(spe,
                          merge,
                          cluster_name = "cluster",
                          new_name = cluster_name,
                          merge_label = NULL,
                          relabel = FALSE,
                          start_from = NULL,
                          sep = "&") {
  if (!cluster_name %in% colnames(SummarizedExperiment::colData(spe)))
    stop("'", cluster_name, "' not found in colData(spe).")
  cl_chr <- as.character(spe[[cluster_name]])

  # Numbering base for relabel. If not given, infer it (e.g. 0 like Seurat, or 1)
  # from the existing numeric labels so relabel preserves what getClusters made.
  if (is.null(start_from)) start_from <- .inferStartFrom(cl_chr)

  # Accept a single vector or a list of vectors.
  if (!is.list(merge)) merge <- list(merge)

  for (grp in merge) {
    grp <- as.character(grp)
    present <- grp[grp %in% cl_chr]
    if (length(present) < 2L) {
      warning("Merge group {", paste(grp, collapse = ", "),
              "} has fewer than two clusters present; skipping.")
      next
    }
    lab <- if (!is.null(merge_label) && length(merge) == 1L) merge_label
           else paste(present, collapse = sep)
    cl_chr[cl_chr %in% present] <- lab
  }

  if (relabel) {
    rl <- .relabelBySize(cl_chr, start_from)
    cl_chr <- rl$labels
    lvls <- rl$levels
  } else {
    lvls <- .orderClusterLevels(cl_chr)
  }

  spe[[new_name]] <- factor(cl_chr, levels = lvls)
  return(spe)
}


#' Renumber all clusters by size.
#'
#' Renumbers the clusters in a colData column \code{1..K} by size (largest
#' first), with "unassigned" kept as the last level. Useful for tidying up the
#' compound labels left by \code{getSubClusters(relabel = FALSE)} /
#' \code{mergeClusters(relabel = FALSE)} after a multi-step workflow.
#'
#' Renumbering is size-based, so run it \strong{before} cell-type annotation -
#' relabelling after a cluster-to-cell-type mapping would scramble that mapping.
#'
#' @param spe A SpatialExperiment object.
#' @param cluster_name Name of the cluster column in
#' \link[SummarizedExperiment]{colData} to renumber.
#' @param new_name Name of the column to store the result. Defaults to
#' overwriting cluster_name.
#' @param start_from Integer at which numbering starts. Defaults to NULL, which
#' infers the base (0- or 1-based) from the existing labels.
#' @return A SpatialExperiment with the renumbered clusters in colData.
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- normalizeAssay(spe)
#' spe <- runPCA(spe)
#' spe <- findNbrsSNN(spe, dimred = "PCA")
#' spe <- getClusters(spe, resolution = 0.5)
#' spe <- getSubClusters(spe, cluster = 2, resolution = 0.5)
#' # Tidy the 2_1/2_2/... labels into plain 1..K numbers by size
#' spe <- relabelClusters(spe)
relabelClusters <- function(spe,
                            cluster_name = "cluster",
                            new_name = cluster_name,
                            start_from = NULL) {
  if (!cluster_name %in% colnames(SummarizedExperiment::colData(spe)))
    stop("'", cluster_name, "' not found in colData(spe).")
  cl_chr <- as.character(spe[[cluster_name]])
  if (is.null(start_from)) start_from <- .inferStartFrom(cl_chr)

  rl <- .relabelBySize(cl_chr, start_from)
  spe[[new_name]] <- factor(rl$labels, levels = rl$levels)
  return(spe)
}


# Infer the numbering base (e.g. 0 like Seurat, or 1) from existing numeric
# cluster labels, so relabelling preserves what getClusters produced.
.inferStartFrom <- function(cl_chr) {
  orig_num <- suppressWarnings(as.numeric(setdiff(unique(cl_chr), "unassigned")))
  if (any(!is.na(orig_num))) min(orig_num, na.rm = TRUE) else 1
}

# Renumber cluster labels by size (largest first), starting at start_from, with
# "unassigned" kept as the last level. Returns the relabelled character vector
# and the corresponding factor levels.
.relabelBySize <- function(cl_chr, start_from) {
  keep <- cl_chr != "unassigned"
  ord  <- names(sort(table(cl_chr[keep]), decreasing = TRUE))
  ids  <- seq_along(ord) - 1L + start_from
  map  <- stats::setNames(ids, ord)
  cl_chr[keep] <- as.character(map[cl_chr[keep]])
  lvls <- c(as.character(ids), if (any(!keep)) "unassigned")
  list(labels = cl_chr, levels = lvls)
}
