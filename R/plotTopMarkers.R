#' Dot plot or heatmap of the top marker genes per cluster.
#'
#' Takes the output of \link[scider]{getMarkers} and visualises the top \code{n}
#' marker genes of each cluster, either as a dot plot or a genes-by-clusters
#' expression heatmap. Genes are block-ordered by cluster (in the order of
#' \code{names(markers)}); a gene that is a top marker for more than one cluster
#' is shown once, under the first cluster it appears in.
#'
#' Two expression profiles are available via \code{profile}:
#' \itemize{
#'   \item \strong{"pseudobulk"} (default): per-cluster log-CPM from
#'     \link[scider]{spe2PB} + \link[edgeR]{cpm} - the same space the markers were
#'     tested in. Dot colour is (scaled) log-CPM; dot size is the proportion of
#'     cells expressing the gene.
#'   \item \strong{"cell"}: Seurat-style cell-level summaries via
#'     \link[scider]{plotDots} (dot plot) or cell-level mean expression (heatmap).
#' }
#'
#' @param markers The named list returned by \link[scider]{getMarkers} (one
#'   data frame per cluster, gene names as row names).
#' @param spe The SpatialExperiment object the markers were computed from.
#' @param n Integer. Number of top markers to take from each cluster. Default 5.
#' @param min.prop Numeric in \[0, 1\]. Drop candidate markers whose detection rate
#'   (proportion of cells expressing) in their cluster is below this value before
#'   taking the top \code{n}, so low-prevalence genes are skipped and back-filled
#'   from further down the ranking. Default 0 (no filtering).
#' @param type Character. "dotplot" (default) or "heatmap".
#' @param profile Character. "pseudobulk" (default) or "cell"; see Details.
#' @param cluster_name Character. Column in \code{colData(spe)} holding the
#'   cluster labels used in \code{markers}. Default "cluster".
#' @param assay Name of the assay to use for expression. Only used when
#'   \code{profile = "cell"} (pseudo-bulk always uses raw counts). Default
#'   "counts".
#' @param scale Logical. Scale expression per gene (z-score across clusters for
#'   the colour). Default TRUE.
#' @param range Numeric. Lower and upper caps for the colour values; values
#'   outside are clamped (winsorised). A single value is taken as the upper
#'   bound. Default NULL picks \code{c(-1.5, 2.2)} when \code{scale = TRUE}
#'   (suited to z-scores) and \code{c(-Inf, Inf)} (no capping) otherwise.
#' @param cols Custom colour palette. For a cell-level dot plot a gradient passed
#'   to \link[scider]{plotDots}; otherwise the full colour vector for the
#'   gradient (dot plot) or \link[pheatmap]{pheatmap} (heatmap).
#' @param dot.scale Numeric. Scales the dot radius. See
#'   \link[ggplot2]{scale_radius}. Default 6.
#' @param ... Further arguments passed to \link[scider]{plotDots} (cell-level
#'   dot plot) or \link[pheatmap]{pheatmap} (heatmap).
#'
#' @return A ggplot object when \code{type = "dotplot"}, or a pheatmap object
#'   when \code{type = "heatmap"}.
#'
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- normalizeAssay(spe)
#' spe <- runPCA(spe)
#' spe <- findNbrsSNN(spe, dimred = "PCA")
#' spe <- getClusters(spe, resolution = 0.5)
#' markers <- getMarkers(spe)
#' plotTopMarkers(markers, spe, n = 5)
#' plotTopMarkers(markers, spe, n = 5, type = "heatmap")
plotTopMarkers <- function(markers,
                           spe,
                           n = 5,
                           min.prop = 0,
                           type = c("dotplot", "heatmap"),
                           profile = c("pseudobulk", "cell"),
                           cluster_name = "cluster",
                           assay = "counts",
                           scale = TRUE,
                           range = NULL,
                           cols = NULL,
                           dot.scale = 6,
                           ...) {
    type <- match.arg(type)
    profile <- match.arg(profile)
    if (is.data.frame(markers) || !is.list(markers) || is.null(names(markers)))
        stop("'markers' must be the named list returned by getMarkers().")
    if (!cluster_name %in% names(SummarizedExperiment::colData(spe)))
        stop("'", cluster_name, "' not found in colData(spe).")
    if (is.null(range)) range <- if (scale) c(-1.5, 2.2) else c(-Inf, Inf)

    # Per-cluster detection rates for the min.prop filter (same quantity as the
    # 'Prop' column from getMarkers); only needed when filtering.
    detprop <- NULL
    if (min.prop > 0) {
        det <- .clusterDetection(spe, cluster_name)
        detprop <- sweep(det$nexpr, 2, det$ncells, "/")
    }

    # Top n markers per cluster, block-ordered by cluster; keep first occurrence.
    # Genes below min.prop in their cluster are dropped before taking the top n,
    # so the slots are back-filled from further down the ranking.
    features <- character(0)
    for (cl in names(markers)) {
        g <- rownames(markers[[cl]])
        if (min.prop > 0 && cl %in% colnames(detprop))
            g <- g[detprop[g, cl] >= min.prop]
        g <- head(g, n)
        features <- c(features, setdiff(g, features))
    }
    features <- features[features %in% rownames(spe)]
    if (length(features) == 0)
        stop("None of the top markers are present in rownames(spe).")

    ## Cell-level profile (Seurat-style): defer to plotDots / cell averaging.
    if (profile == "cell") {
        if (type == "dotplot")
            return(plotDots(spe, feature = features, assay = assay,
                            group.by = cluster_name, scale = scale,
                            range = range,
                            cols = cols, dot.scale = dot.scale, ...))
        exprs <- as.matrix(
            SummarizedExperiment::assay(spe, assay)[features, , drop = FALSE])
        group <- factor(spe[[cluster_name]])
        clusters <- names(markers)[names(markers) %in% levels(group)]
        mat <- vapply(clusters, function(g)
            log1p(rowMeans(expm1(exprs[, group == g, drop = FALSE]))),
            numeric(length(features)))
        rownames(mat) <- features
        if (scale) { mat <- t(scale(t(mat))); mat[is.na(mat)] <- 0 }
        mat <- .capLimits(mat, range)
        return(.markerHeatmap(mat, cols = cols, scaled = scale, ...))
    }

    ## Pseudo-bulk profile: per-cluster log-CPM (same space as getMarkers()).
    y <- spe2PB(spe, group.id = cluster_name)
    clusters <- names(markers)[names(markers) %in% colnames(y)]
    mat <- edgeR::cpm(y, log = TRUE)[features, clusters, drop = FALSE]
    if (scale) { mat <- t(scale(t(mat))); mat[is.na(mat)] <- 0 }
    mat <- .capLimits(mat, range)

    if (type == "heatmap")
        return(.markerHeatmap(mat, cols = cols, scaled = scale, ...))

    # Dot plot: colour = (scaled) log-CPM, size = proportion of cells expressing.
    det <- .clusterDetection(spe, cluster_name)
    prop <- sweep(det$nexpr, 2, det$ncells, "/")[features, clusters, drop = FALSE]
    long <- data.frame(
        cluster    = factor(rep(clusters, each = length(features)),
                            levels = clusters),
        feature    = factor(rep(features, times = length(clusters)),
                            levels = features),
        average    = as.vector(mat),
        percentage = as.vector(prop))
    p <- ggplot(long) +
        geom_point(aes(.data[["cluster"]], .data[["feature"]],
                       size = .data[["percentage"]],
                       fill = .data[["average"]]),
                   shape = 21, colour = "black", stroke = 0.3) +
        theme_minimal() +
        theme(axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_radius(range = c(0, dot.scale)) +
        labs(x = cluster_name, y = NULL, size = "Prop.",
             fill = if (scale) "z-score" else "logCPM")
    if (is.null(cols)) p <- p + scale_fill_gradient(low = "lightgrey",
                                                    high = "darkred")
    else p <- p + scale_fill_gradientn(colours = cols)
    p
}

# Winsorise a matrix to the given colour range (length-2; a single value is
# taken as the upper bound). Inf/-Inf leave the data unchanged.
.capLimits <- function(mat, lim) {
    if (length(lim) == 1) lim <- c(-Inf, lim)
    mat[mat < lim[1]] <- lim[1]
    mat[mat > lim[2]] <- lim[2]
    mat
}

# pheatmap of a genes x clusters matrix, with the marker block order preserved
# (no row/column clustering) and a symmetric diverging scale when row-scaled.
.markerHeatmap <- function(mat, cols = NULL, scaled = TRUE, ...) {
    if (is.null(cols))
        cols <- rev(grDevices::colorRampPalette(c(
            "#ED254EFF", "#EF6079FF", "#F1F4FFFF",
            "#97B3D0FF", "#011936FF"))(100))
    breaks <- NA
    if (scaled) {
        mx <- max(abs(mat))
        breaks <- seq(-mx, mx, length.out = length(cols) + 1)
    }
    pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
                       color = cols, breaks = breaks, border_color = "white",
                       angle_col = 45, ...)
}
