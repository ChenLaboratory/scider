#' Find up-regulated marker genes for all clusters (quasi-NB GLM).
#'
#' Pseudo-bulks cells by cluster via \link[scider]{spe2PB} and, for each cluster,
#' tests its genes against the average of all other clusters (1-vs-rest) using a
#' quasi-likelihood approach. Because the full design is saturated (one
#' pseudo-bulk per cluster, 0 residual df), each contrast is tested by fitting
#' the corresponding null design (one fewer column, 1 residual df) with
#' \link[edgeR]{glmQLFit} and using the quasi-deviance chi-square statistic.
#' The test is \strong{one-sided}: p-values reflect up-regulation in the cluster,
#' so each table ranks that cluster's positive markers.
#'
#' Each per-cluster table holds \strong{all} genes, sorted by \code{sort.by} with
#' an adjusted p-value column added via \code{adjust.method}. Use
#' \link[scider]{topMarkers} to extract and combine the top genes per cluster.
#' For a direct comparison between two clusters (or two groups of clusters), see
#' \link[scider]{getDE}.
#'
#' @param spe A SpatialExperiment object.
#' @param cluster_name Character. Column name in \code{colData(spe)} containing
#'   cluster labels. Default \code{"cluster"}.
#' @param dispersion Numeric. Fixed NB dispersion (BCV^2) used for the GLM fits.
#'   Required because the saturated full design (one pseudo-bulk per cluster)
#'   leaves no residual df for dispersion estimation. Default 0.05
#'   (BCV \eqn{\approx} 0.224).
#' @param method Character. Testing pipeline: "quasi" (default) fits each null
#'   design with \link[edgeR]{glmQLFit} and uses the quasi-deviance chi-square
#'   statistic; "lrt" uses the standard \link[edgeR]{glmFit} +
#'   \link[edgeR]{glmLRT} likelihood-ratio test. Both use the preset
#'   \code{dispersion}.
#' @param adjust.method Character. Multiple-testing adjustment method passed to
#'   \link[stats]{p.adjust} (e.g. "BH", "BY", "holm", "bonferroni", "none").
#'   Default "BH".
#' @param sort.by Character. How to order each table: "PValue" (default),
#'   "logFC" (by absolute log fold change), "stat", or "none".
#' @param ... Additional arguments passed to \link[edgeR]{glmQLFit} for the
#'   null-model fits (used only when \code{method = "quasi"}).
#'
#' @return A named \code{list} of data frames, one per cluster, each containing
#'   all genes. Columns:
#'   \describe{
#'     \item{logFC}{log2 fold change of the cluster vs the average of the rest.}
#'     \item{logCPM}{average log2 counts-per-million across pseudo-bulk samples.}
#'     \item{Prop}{proportion of cells in the cluster with a non-zero count for
#'       the gene.}
#'     \item{stat}{signed quasi statistic (signed square root of the quasi
#'       chi-square statistic); positive means up-regulated in the cluster.}
#'     \item{PValue}{one-sided (up-regulation) p-value.}
#'     \item{FDR}{adjusted p-value (column named per \code{adjust.method}; "FDR"
#'       for BH/BY, "FWER" for family-wise methods). Absent when
#'       \code{adjust.method = "none"}.}
#'   }
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
getMarkers <- function(spe,
                       cluster_name = "cluster",
                       dispersion = 0.05,
                       method = c("quasi", "lrt"),
                       adjust.method = "BH",
                       sort.by = "PValue",
                       ...) {
    method <- match.arg(method)
    pb <- .pbClusterFit(spe, cluster_name, dispersion)

    # 1-vs-rest contrast matrix (column i: cluster i vs average of the rest)
    contr <- matrix(-1 / (pb$n - 1), pb$n, pb$n)
    diag(contr) <- 1

    res <- lapply(seq_len(pb$n), function(i) {
        tab <- .contrastTest(pb, contr[, i], method = method,
                             alternative = "up", dispersion = dispersion, ...)
        tab$Prop <- pb$nexpr[rownames(tab), i] / pb$ncells[i]
        tab <- tab[, c("logFC", "logCPM", "Prop", "stat", "PValue")]
        # All genes, sorted with an adjusted-p column; selection is topMarkers'.
        .topTags(tab, n = Inf, adjust.method = adjust.method,
                 sort.by = sort.by, p.value = 1)
    })
    names(res) <- pb$clst_names
    res
}


#' Combine the top markers of each cluster into one data frame.
#'
#' Convenience wrapper around \link[scider]{getMarkers} output: takes the top
#' \code{n} genes from each cluster's (already ranked) table and stacks them into
#' a single data frame, with a leading \code{cluster} column and a \code{gene}
#' column (the same gene can be a top marker for more than one cluster, so row
#' names cannot stay unique). Analogous to the single table returned by
#' \code{Seurat::FindAllMarkers}.
#'
#' @param markers The named list returned by \link[scider]{getMarkers}.
#' @param n Integer. Number of top genes to take from each cluster. Default 10.
#' @param p.value Numeric. Keep only genes whose adjusted p-value (or raw
#'   \code{PValue} if no adjusted column is present) is at or below this cutoff,
#'   applied before taking the top \code{n}. Default 1 (no filtering).
#' @param min.prop Numeric in \[0, 1\]. Drop genes whose detection rate (the
#'   \code{Prop} column - proportion of cells expressing in the cluster) is below
#'   this value before taking the top \code{n}, so the slots are back-filled from
#'   further down the ranking. Default 0 (no filtering).
#' @return A data frame with columns \code{cluster}, \code{gene}, and the
#'   per-cluster statistics from \link[scider]{getMarkers}, ordered by cluster
#'   then by each cluster's ranking.
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- normalizeAssay(spe)
#' spe <- runPCA(spe)
#' spe <- findNbrsSNN(spe, dimred = "PCA")
#' spe <- getClusters(spe, resolution = 0.5)
#' markers <- getMarkers(spe)
#' top <- topMarkers(markers, n = 10)
topMarkers <- function(markers, n = 10, p.value = 1, min.prop = 0) {
    if (is.data.frame(markers) || !is.list(markers) || is.null(names(markers)))
        stop("'markers' must be the named list returned by getMarkers().")

    out <- lapply(names(markers), function(cl) {
        tab <- markers[[cl]]
        if (p.value < 1) {
            adj <- intersect(c("FDR", "FWER"), colnames(tab))
            sig <- if (length(adj)) tab[[adj[1]]] else tab$PValue
            tab <- tab[sig <= p.value, , drop = FALSE]
        }
        if (min.prop > 0 && "Prop" %in% colnames(tab))
            tab <- tab[tab$Prop >= min.prop, , drop = FALSE]
        tab <- head(tab, n)
        if (nrow(tab) == 0) return(NULL)
        data.frame(cluster = cl, gene = rownames(tab), tab,
                   row.names = NULL, check.names = FALSE)
    })
    res <- do.call(rbind, out)
    if (!is.null(res))
        res$cluster <- factor(res$cluster, levels = names(markers))
    res
}


#' Differential expression between two clusters or two groups of clusters.
#'
#' Pseudo-bulks cells by cluster via \link[scider]{spe2PB} and performs a
#' \strong{two-sided} quasi-likelihood test comparing the average of one set of
#' clusters (\code{cluster}) against another (\code{cluster_2}), analogous to
#' \link[edgeR]{glmQLFTest}. The contrast is the difference of the two group
#' means, tested by fitting the null design with \link[edgeR]{glmQLFit} (see
#' \link[scider]{getMarkers} for why the null design is used). The result table
#' is post-processed like \link[edgeR]{topTags}.
#'
#' @param spe A SpatialExperiment object.
#' @param cluster Character vector. One or more cluster labels forming the first
#'   group. Positive \code{logFC} means up-regulated in this group.
#' @param cluster_2 Character vector. One or more cluster labels forming the
#'   second (reference) group. Must not overlap \code{cluster}.
#' @param cluster_name Character. Column name in \code{colData(spe)} containing
#'   cluster labels. Default \code{"cluster"}.
#' @param dispersion Numeric. Fixed NB dispersion (BCV^2) used for the GLM fits.
#'   Default 0.05 (BCV \eqn{\approx} 0.224).
#' @param method Character. Testing pipeline: "quasi" (default,
#'   \link[edgeR]{glmQLFit} quasi-deviance test) or "lrt" (\link[edgeR]{glmFit} +
#'   \link[edgeR]{glmLRT}). Both use the preset \code{dispersion}.
#' @param n Integer. Maximum number of genes to return. Default \code{Inf}
#'   (all genes).
#' @param adjust.method Character. Multiple-testing adjustment method passed to
#'   \link[stats]{p.adjust}. Default "BH".
#' @param sort.by Character. How to order the table: "PValue" (default),
#'   "logFC", "stat", or "none".
#' @param p.value Numeric. Cutoff on the adjusted p-value; genes above it are
#'   dropped. Default 1 (no filtering).
#' @param ... Additional arguments passed to \link[edgeR]{glmQLFit} for the
#'   null-model fit (used only when \code{method = "quasi"}).
#'
#' @return A data frame with columns \code{logFC} (log2 fold change of
#'   \code{cluster} vs \code{cluster_2}), \code{logCPM}, \code{Prop.1} and
#'   \code{Prop.2} (proportion of cells expressing the gene in group 1 and
#'   group 2 respectively), \code{stat} (signed quasi statistic), \code{PValue}
#'   (two-sided), and the adjusted p-value column.
#'
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' spe <- normalizeAssay(spe)
#' spe <- runPCA(spe)
#' spe <- findNbrsSNN(spe, dimred = "PCA")
#' spe <- getClusters(spe, resolution = 0.5)
#' # cluster 1 vs cluster 2
#' de <- getDE(spe, cluster = 1, cluster_2 = 2)
#' # clusters 1 & 3 vs clusters 2 & 4
#' de2 <- getDE(spe, cluster = c(1, 3), cluster_2 = c(2, 4))
getDE <- function(spe,
                   cluster,
                   cluster_2,
                   cluster_name = "cluster",
                   dispersion = 0.05,
                   method = c("quasi", "lrt"),
                   n = Inf,
                   adjust.method = "BH",
                   sort.by = "PValue",
                   p.value = 1,
                   ...) {
    method <- match.arg(method)
    pb <- .pbClusterFit(spe, cluster_name, dispersion)

    g1 <- as.character(cluster)
    g2 <- as.character(cluster_2)
    bad <- setdiff(c(g1, g2), pb$clst_names)
    if (length(bad))
        stop("cluster label(s) not found: ", paste(bad, collapse = ", "))
    if (length(intersect(g1, g2)))
        stop("'cluster' and 'cluster_2' must not share any clusters.")

    # Contrast: mean of group 1 minus mean of group 2 (each group averaged).
    cv <- numeric(pb$n)
    cv[match(g1, pb$clst_names)] <-  1 / length(g1)
    cv[match(g2, pb$clst_names)] <- -1 / length(g2)

    tab <- .contrastTest(pb, cv, method = method, alternative = "two.sided",
                         dispersion = dispersion, ...)
    tab$Prop.1 <- rowSums(pb$nexpr[rownames(tab), g1, drop = FALSE]) /
        sum(pb$ncells[g1])
    tab$Prop.2 <- rowSums(pb$nexpr[rownames(tab), g2, drop = FALSE]) /
        sum(pb$ncells[g2])
    tab <- tab[, c("logFC", "logCPM", "Prop.1", "Prop.2", "stat", "PValue")]
    .topTags(tab, n = n, adjust.method = adjust.method,
             sort.by = sort.by, p.value = p.value)
}


# ---- internal helpers -------------------------------------------------------

# Pseudo-bulk by cluster and fit the saturated full GLM. The full design has one
# coefficient per cluster (0 residual df), so it is fitted with glmFit purely to
# recover per-cluster coefficients for logFC; the quasi-dispersion squeeze is
# done per contrast on the null design in .contrastQLTest().
.pbClusterFit <- function(spe, cluster_name, dispersion) {
    if (!requireNamespace("edgeR", quietly = TRUE))
        stop("edgeR is required but is not installed (or can't be loaded)")
    if (!cluster_name %in% names(SummarizedExperiment::colData(spe)))
        stop("'", cluster_name, "' not found in colData(spe).")

    y <- spe2PB(spe, group.id = cluster_name)
    n <- ncol(y)
    clst_names <- colnames(y)
    design <- stats::model.matrix(~ 0 + factor(seq_len(n)))
    fit_full <- edgeR::glmFit(y, design = design, dispersion = dispersion)

    # Per-cluster detection rates, aligned to y's columns.
    det <- .clusterDetection(spe, cluster_name)
    list(y = y, design = design, fit_full = fit_full,
         clst_names = clst_names, n = n,
         logCPM = edgeR::aveLogCPM(y),
         nexpr = det$nexpr[, clst_names, drop = FALSE],
         ncells = det$ncells[clst_names])
}

# For each gene, the number of cells with a non-zero count in each cluster
# (genes x clusters) and the number of cells per cluster. Used for the
# proportion-expressing statistics in getMarkers()/getDE() and plotTopMarkers().
.clusterDetection <- function(spe, cluster_name) {
    if (!cluster_name %in% names(SummarizedExperiment::colData(spe)))
        stop("'", cluster_name, "' not found in colData(spe).")
    counts <- spe@assays@data$counts
    cl <- factor(as.character(SummarizedExperiment::colData(spe)[[cluster_name]]))
    gmat <- stats::model.matrix(~ 0 + cl)
    colnames(gmat) <- gsub("^cl", "", colnames(gmat))
    nexpr <- as.matrix(((counts > 0) * 1) %*% gmat)
    list(nexpr = nexpr, ncells = colSums(gmat))
}

# Test a single 1-df contrast cv at the preset dispersion. Two pipelines:
#   "quasi" - fit the null design (full design with the cv direction removed,
#     n-1 columns -> 1 residual df) with glmQLFit so squeezeVar has df, then use
#     the quasi-deviance chi-square statistic.
#   "lrt"   - standard glmFit + glmLRT likelihood-ratio test on the full fit.
# In both, the statistic is ~chi-square on 1 df, so stat = sign(logFC)*sqrt(chi).
# alternative = "up" gives a one-sided (up-regulation) p-value; "two.sided" the
# usual two-sided p-value. Returns an unsorted, unadjusted table (logFC, logCPM,
# stat, PValue); .topTags() does the rest.
.contrastTest <- function(pb, cv, method = c("quasi", "lrt"),
                          alternative = c("two.sided", "up"),
                          dispersion, ...) {
    method <- match.arg(method)
    alternative <- match.arg(alternative)

    logFC <- drop(pb$fit_full$coefficients %*% cv) / log(2)

    if (method == "quasi") {
        # Basis for the null space of the contrast (drop the cv direction).
        Q <- qr.Q(qr(matrix(cv, ncol = 1L)), complete = TRUE)
        design0 <- pb$design %*% Q[, -1L, drop = FALSE]
        fit1 <- edgeR::glmQLFit(pb$y, design0, dispersion = dispersion, ...)
        chi <- fit1$deviance.adj / fit1$average.ql.dispersion
        df  <- fit1$df.residual.adj
    } else {
        lrt <- edgeR::glmLRT(pb$fit_full, contrast = cv)
        chi <- lrt$table$LR
        df  <- lrt$df.test
    }

    stat  <- sign(logFC) * sqrt(chi)
    p_two <- stats::pchisq(chi, df = df, lower.tail = FALSE)
    pval  <- if (alternative == "up")
        ifelse(logFC >= 0, p_two / 2, 1 - p_two / 2) else p_two

    data.frame(logFC  = logFC,
               logCPM = pb$logCPM,
               stat   = stat,
               PValue = pval,
               row.names = rownames(pb$y))
}

# Post-process a raw test table like edgeR::topTags: add the adjusted p-value
# column, order by sort.by, drop genes failing the p.value cutoff, then keep the
# top n. The FDR column is computed on the full gene set before any filtering.
.topTags <- function(tab, n = 10, adjust.method = "BH",
                     sort.by = "PValue", p.value = 1) {
    sort.by <- match.arg(sort.by, c("PValue", "logFC", "stat", "none"))

    # Adjusted p-value column (named as edgeR does).
    adj.name <- NULL
    if (adjust.method != "none") {
        adj.name <- if (adjust.method %in% c("BH", "fdr", "BY")) "FDR"
            else if (adjust.method %in%
                     c("holm", "hochberg", "hommel", "bonferroni")) "FWER"
            else adjust.method
        tab[[adj.name]] <- stats::p.adjust(tab$PValue, method = adjust.method)
    }

    # Order. Break PValue ties by descending |stat|, since many top genes can
    # share PValue = 0 (underflow) and would otherwise be in arbitrary order.
    o <- switch(sort.by,
        PValue = order(tab$PValue, -abs(tab$stat)),
        logFC  = order(abs(tab$logFC), decreasing = TRUE),
        stat   = order(abs(tab$stat),  decreasing = TRUE),
        none   = seq_len(nrow(tab)))
    tab <- tab[o, , drop = FALSE]

    # Filter by significance cutoff.
    if (p.value < 1) {
        sig <- if (is.null(adj.name)) tab$PValue else tab[[adj.name]]
        tab <- tab[sig <= p.value, , drop = FALSE]
    }

    # Keep top n.
    if (n < nrow(tab)) tab <- tab[seq_len(n), , drop = FALSE]
    tab
}
