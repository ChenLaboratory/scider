#' Test for density correlation between two cell types.
#'
#' @param spe A SpatialExperiment object.
#' @param coi Character vector for cell types of interest for density
#' correlation analysis. Default is NULL, which is to consider all cell types
#' previously calculated in the gridDensity() step.
#' @param roi Character. The name of the group or cell type on which
#' the roi is computed. Default is NULL for no subsetting cell types by ROI
#' @param probs A numeric scalar. The threshold of proportion that used to
#' filter grids by density when ROIs have not been identified previously.
#' Ignored if 'roi' is present in the 'metadata' component of spe. Default to 0.85.
#' @param min.density A numeric value. The cut-off on the same density score
#' as 'probs': grids with a score of at least 'min.density' are kept. Ignored if 
#' 'roi' is given. Default is NULL and overwrites probs.
#' @param trace Logical. If TRUE, print the process of testing. Default to FALSE.
#' @param method Character. The testing method, either "permutation" (default)
#' or "modified.ttest". See Details.
#' @param id Character. The column of colData(spe) holding the cell type
#' labels that were used in gridDensity(). Only used by method "permutation".
#' Default to "cell_type".
#' @param nperm Integer. The number of random relabellings used by method
#' "permutation". Default to 999.
#' @param permute.within Character or NULL. The column of colData(spe) that
#' splits cells into groups (e.g. tissue compartments) for method
#' "permutation". Cell type labels are then shuffled only among cells of the
#' same group. Default to NULL, which shuffles labels across all cells. See
#' Details.
#' @param cpu_threads Integer. The number of threads used by method
#' "permutation". The results do not depend on it. Default to 1.
#'
#' @details
#' Method "permutation" tests each pair of cell types against random
#' labelling. Every cell keeps its location and the cell type labels are
#' shuffled across all cells of the sample, after which the cell type
#' densities are recomputed with the same kernel and bandwidth as in
#' gridDensity(). The statistic is the Spearman correlation between the two
#' densities across grids: within each ROI for the ROI level results, and
#' pooled over all ROIs after centering the ranks at the ROI level for the
#' overall results. The column 'cor.expected' gives the mean of the
#' statistic under random labelling, which is the correlation expected from
#' the shared dependence of both densities on local cellularity alone.
#' 'p.Pos' and 'p.Neg' are one sided Monte Carlo p-values, (1 + b) / (1 +
#' nperm), for the observed correlation being larger or smaller than
#' expected under random labelling. The 'shift' (cor.coef - cor.expected) 
#' measures the strength of co-localization (positive) or
#' segregation (negative) beyond random labelling. The overall results also
#' give a per pair summary of the ROI level results: 'n.roi' ROIs
#' with a defined correlation, 'n.above' and 'n.below' ROIs with the observed
#' correlation above or below its expected value, and 'n.sig.Pos' and
#' 'n.sig.Neg' ROIs with p.Pos or p.Neg below 0.05. With many grids, small
#' shifts reach significance, so rank pairs by the shift and check the ROI
#' summary for consistency. Use set.seed() for reproducible p-values; 
#' the results do not depend on 'cpu_threads'.
#'
#' With 'permute.within', the null hypothesis is random labelling within
#' groups: each group keeps its cell locations and cell type counts, and only
#' the labels inside a group are shuffled. Structure already encoded by the
#' groups, such as the split between tumour and stroma compartments, then
#' becomes part of the expected value, and the test asks whether cell types
#' are organized beyond it. For two cell types in the same group a positive
#' shift means they share neighbourhoods more than two random cells of that
#' group would; for cell types in different groups it means they sit closer
#' to each other than random cells of their groups would. Groups should be
#' coarse: a cell type that is the only type in every group it occupies is
#' never moved, so pairs involving it are only tested through the other cell
#' type.
#'
#' Method "modified.ttest" is the original approach: a Pearson correlation
#' within each ROI tested by the modified t-test of Dutilleul (1993) from
#' SpatialPack, with the ROI level p-values combined by Fisher's method. Its
#' effective degrees of freedom are estimated from the spatial correlogram
#' within each ROI.
#'
#' @return A list of two DataFrames, 'ROI' (results per ROI) and 'overall'
#' (results combined across ROIs). When 'roi' is NULL the selected grids
#' form a single ROI, labelled "1", and the overall results equal the ROI
#' level results. For method
#' "permutation" the columns are celltype1, celltype2, ROI, ngrid, cor.coef
#' (observed Spearman correlation), cor.expected (mean under random
#' labelling), shift, p.Pos and p.Neg; the overall DataFrame has n.roi,
#' n.above, n.below, n.sig.Pos and n.sig.Neg (see Details) in place of ROI
#' and ngrid. For method "modified.ttest" the columns are
#' celltype1, celltype2, ROI, ngrid, cor.coef, t, df, p.Pos and p.Neg.
#' @export
#'
#' @importFrom SpatialPack modified.ttest
#' @importFrom stats pt
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
#'
#' spe <- gridDensity(spe, coi = coi)
#'
#' spe <- findROI(spe, coi = coi, method = "walktrap")
#'
#' set.seed(1)
#' result <- corDensity(spe, roi = coi, nperm = 199)
#'
corDensity <- function(spe, coi = NULL, roi = NULL, probs = 0.85,
                       min.density = NULL, trace = FALSE,
                       method = c("permutation", "modified.ttest"),
                       id = "cell_type", nperm = 999,
                       permute.within = NULL, cpu_threads = 1) {
  method <- match.arg(method)
  if (!is.null(permute.within) && method != "permutation") {
    warning("'permute.within' is only used by method = \"permutation\" ",
            "and is ignored.")
  }
  if (!("grid_density" %in% names(spe@metadata))) {
    stop("Please run gridDensity before using this function.")
  }

  dens_dat <- as.data.frame(spe@metadata$grid_density)

  # get cell type info
  den_cols <- colnames(dens_dat)[grepl("density_", colnames(dens_dat))]
  den_cols <- den_cols[den_cols != "density_overall"]

  if(any(!is.null(coi))){
    coi_clean <- paste0("density_", janitor::make_clean_names(coi))
    coi.exist <- coi_clean %in% den_cols
    if(any(!coi.exist))
      stop(paste0(paste(coi[!coi.exist], collapse=", "), " not found in the data."))
    den_cols <- den_cols[den_cols %in% coi_clean]
  }

  nCT <- length(den_cols)
  if (nCT < 2) stop("Please run gridDensity for at least two of the cell types specified in 'coi'.")

  # construct data table
  if (is.null(roi)) {
    if (method == "permutation") {
      if (!("density_overall" %in% colnames(dens_dat))) {
        stop("density_overall is not found. Please rerun gridDensity().")
      }
      sel_score <- dens_dat$density_overall
    } else {
      sel_score <- rowMeans(as.matrix(dens_dat[, which(colnames(dens_dat) %in% den_cols), drop = FALSE]))
    }
    if (!is.null(min.density)) {
      if (!is.numeric(min.density) || length(min.density) != 1L ||
          is.na(min.density)) {
        stop("'min.density' must be a single number.")
      }
      if (!missing(probs)) {
        message("Overwriting the probs argument. Grids are filtered by the ",
                "min.density value.")
      }
      kp <- sel_score >= min.density
      if (!any(kp)) {
        stop("No grid has a density of at least 'min.density' (the largest ",
             "is ", signif(max(sel_score), 3), ").")
      }
    } else {
      kp <- sel_score >= quantile(sel_score, probs = probs)
    }
    dens_dat_filter <- dens_dat[kp, ]
    rois <- data.frame(component=gl(1,sum(kp)),
                       members=dens_dat_filter$node,
                       x=dens_dat_filter$node_x,
                       y=dens_dat_filter$node_y)
  } else {
    if (!is.null(min.density)) {
      warning("'min.density' is only used when 'roi' is NULL and is ignored.")
    }
    roi <- gsub("_roi$", "", roi)
    roi <- janitor::make_clean_names(roi)
    roi <- paste(c(sort(roi),"roi"), collapse="_")
    if (is.null(spe@metadata[[roi]])) {
      stop(paste(
        roi, " is not found in metadata of spe. Please run findROI() first."
      ))
    }
    rois <- as.data.frame(spe@metadata[[roi]])
  }

  model_data <- merge(rois, dens_dat,
                      by.x = "members",
                      by.y = "node", all.x = TRUE, sort = FALSE
  )

  if (method == "modified.ttest") {
    res <- .corDensityModTtest(model_data, den_cols, trace)
  } else {
    res <- .corDensityPerm(spe, dens_dat, model_data, den_cols, id, nperm,
                           trace, permute.within, cpu_threads)
  }

  result.ROI <- S4Vectors::DataFrame(res$ROI)
  result.overall <- S4Vectors::DataFrame(res$overall)

  list(ROI = result.ROI, overall = result.overall)
}

.gridLattice <- function(x, y) {
  tol <- 1e-6 * max(diff(range(x)), diff(range(y)), 1)
  # smallest spacing of u among nodes that share v
  minStep <- function(u, v) {
    vk <- round(v / tol)
    o <- order(vk, u)
    d <- diff(u[o])
    d <- d[diff(vk[o]) == 0 & d > tol]
    if (length(d)) min(d) else NA_real_
  }
  gcd <- function(v) {
    g <- 0
    for (z in unique(v)) {
      while (z > 0) { t <- g %% z; g <- z; z <- t }
    }
    max(g, 1)
  }
  toIndex <- function(u, step) {
    if (is.na(step)) return(list(i = integer(length(u)), s = 1))
    i <- round((u - min(u)) / step)
    if (max(abs(min(u) + i * step - u)) > 1e-3 * step) {
      stop("The grid nodes are not on a regular lattice. ",
           "Please rerun gridDensity().")
    }
    g <- gcd(i)
    list(i = as.integer(i / g), s = step * g)
  }
  ax <- toIndex(x, minStep(x, y) / 2)
  ry <- toIndex(y, minStep(y, x) / 2)
  if (is.na(minStep(y, x))) ry$s <- ax$s
  if (is.na(minStep(x, y))) ax$s <- ry$s
  par <- unique((ax$i + ry$i) %% 2L)
  hex <- length(par) == 1L && length(unique(ax$i %% 2L)) == 2L
  list(a = ax$i, r = ry$i, sa = ax$s, sr = ry$s,
       cs = if (hex) 2L else 1L, parity = if (hex) par else 0L)
}

# Cell type name for the output, e.g. "density_b_cells" -> "B cells".
.ctName <- function(ct) {
  janitor::make_clean_names(ct, case = "sentence", replace = c("density" = ""))
}

# Original method: Pearson correlation per ROI with the modified t-test of
# SpatialPack, combined across ROIs by Fisher's method.
.corDensityModTtest <- function(model_data, den_cols, trace) {
  nCT <- length(den_cols)
  component <- factor(model_data$component)
  nGrids <- table(component)
  cpnts <- names(nGrids)
  nCpnts <- length(cpnts)

  nRows <- choose(nCT, 2) * nCpnts
  result.ROI <- data.frame(
    "celltype1" = rep("", nRows),
    "celltype2" = rep("", nRows),
    "ROI" = 0,
    "ngrid" = 0,
    "cor.coef" = 0,
    "t" = 0,
    "df" = 0,
    "p.Pos" = 0,
    "p.Neg" = 0
  )

  result.overall <- result.ROI[seq_len(choose(nCT, 2)), c(1, 2, 5, 8, 9)]

  for (i in seq_len((nCT - 1))) {
    for (j in (i + 1):nCT) {
      ct1 <- den_cols[i]
      ct2 <- den_cols[j]
      n1 <- .ctName(ct1)
      n2 <- .ctName(ct2)
      m <- choose(nCT, 2) - (choose(nCT - j, 1) + choose(nCT - i, 2))

      for (k in seq_len(nCpnts)) {
        if (trace) cat(paste("i =", i, ", j =", j, ", ROI", k, "\n"))
        data <- model_data[model_data$component == cpnts[k], c(ct1, ct2, "x", "y")]

        res <- modified.ttest(x=data[,1], y=data[,2],
                              coords=data[,c("x","y")], nclass=7)
        if(res$dof < 0) {
          res <- modified.ttest(x=data[,1], y=data[,2],
                                coords=data[,c("x","y")], nclass=1)
          res$dof <- 1
        }
        # SpatialPack returns Fstat = r^2 / (1 - r^2), so this is the usual t
        tstat <- sqrt( res$Fstat * res$dof ) * sign(res$corr)

        ind <- (m - 1) * nCpnts + k
        result.ROI[ind, "celltype1"] <- n1
        result.ROI[ind, "celltype2"] <- n2
        result.ROI[ind, "ROI"] <- cpnts[k]
        result.ROI[ind, "ngrid"] <- nGrids[k]
        result.ROI[ind, "cor.coef"] <- res$corr
        result.ROI[ind, "t"] <- tstat
        result.ROI[ind, "df"] <- res$dof
        result.ROI[ind, "p.Pos"] <- pt(tstat, df=res$dof, lower.tail = FALSE)
        result.ROI[ind, "p.Neg"] <- pt(tstat, df=res$dof, lower.tail = TRUE)
      }

      # summarize test results across ROIs
      sel <- (m - 1) * nCpnts + seq_len(nCpnts)
      out.cor.coef <- result.ROI$cor.coef[sel]
      out.p.Pos <- result.ROI$p.Pos[sel]
      out.p.Neg <- result.ROI$p.Neg[sel]

      result.overall[m, "celltype1"] <- n1
      result.overall[m, "celltype2"] <- n2
      result.overall[m, "cor.coef"] <- tanh(mean(atanh(out.cor.coef)))
      result.overall[m, "p.Pos"] <- pchisq(-2 * sum(log(out.p.Pos)),
                                           df = 2 * nCpnts, lower.tail = FALSE)
      result.overall[m, "p.Neg"] <- pchisq(-2 * sum(log(out.p.Neg)),
                                           df = 2 * nCpnts, lower.tail = FALSE)
    }
  }
  list(ROI = result.ROI, overall = result.overall)
}

# Random labelling permutation test. The relabelling, densities, ranks and
# correlations of every permutation run in C++ (src/C_corDensityPerm.cpp).
.corDensityPerm <- function(spe, dens_dat, model_data, den_cols, id, nperm,
                            trace, permute.within = NULL, cpu_threads = 1L) {
  nperm <- suppressWarnings(as.integer(nperm))
  if (length(nperm) != 1L || is.na(nperm) || nperm < 1L) {
    stop("'nperm' must be a positive integer.")
  }
  gi <- spe@metadata$grid_info
  if (!is.null(gi$kernel) && !identical(gi$kernel, "gaussian")) {
    stop("method = \"permutation\" only supports the gaussian kernel. ",
         "Please use method = \"modified.ttest\".")
  }
  bw <- gi$bandwidth
  if (is.null(bw)) stop("Bandwidth not found. Please rerun gridDensity().")
  if (is.null(id) || !(id %in% colnames(SummarizedExperiment::colData(spe)))) {
    stop("Column '", id, "' is not found in colData(spe). ",
         "method = \"permutation\" needs the cell type labels used in ",
         "gridDensity(); set 'id' or use method = \"modified.ttest\".")
  }

  # match density columns to cell type labels
  lab <- as.character(spe[[id]])
  lab[is.na(lab)] <- ".NA"
  lv <- unique(lab)
  lv_col <- paste0("density_", vapply(lv, janitor::make_clean_names, character(1)))
  ct_lab <- lv[match(den_cols, lv_col)]
  if (anyNA(ct_lab)) {
    stop(paste(den_cols[is.na(ct_lab)], collapse = ", "),
         " cannot be matched to the labels in colData(spe)$", id, ". ",
         "method = \"permutation\" requires cell type densities; ",
         "use method = \"modified.ttest\" for feature densities.")
  }
  nCT <- length(den_cols)
  ct_idx <- match(lab, ct_lab)                      # NA for other cell types
  ncell <- length(lab)

  # labels are shuffled within groups (a single group by default)
  if (is.null(permute.within)) {
    gid <- rep(1L, ncell)
  } else {
    if (length(permute.within) != 1L ||
        !(permute.within %in% colnames(SummarizedExperiment::colData(spe)))) {
      stop("'permute.within' must be the name of a column of colData(spe).")
    }
    st <- spe[[permute.within]]
    if (anyNA(st)) stop("colData(spe)$", permute.within, " must not contain NA.")
    st <- as.character(st)
    st_idx <- split(seq_len(ncell), st)
    fixed <- ct_lab[vapply(seq_along(ct_lab), function(k) {
      ss <- unique(st[lab == ct_lab[k]])
      all(vapply(ss, function(z) all(lab[st_idx[[z]]] == ct_lab[k]), logical(1)))
    }, logical(1))]
    if (length(fixed)) {
      warning(paste(fixed, collapse = ", "), " is the only cell type in every ",
              "group it occupies and is never permuted.")
    }
    gid <- as.integer(factor(st))
  }

  model_data <- model_data[!is.na(model_data$x_grid), , drop = FALSE]
  comp <- factor(model_data$component)
  cpnts <- levels(comp)
  nCpnts <- length(cpnts)
  comp_int <- as.integer(comp)
  comp_n <- tabulate(comp_int, nbins = nCpnts)

  # lattice coordinates of the grid nodes; cells are counted at their
  # nearest node
  lat <- .gridLattice(dens_dat$x_grid, dens_dat$y_grid)
  xy_node <- cbind(dens_dat$x_grid, dens_dat$y_grid)
  xy_cell <- as.matrix(SpatialExperiment::spatialCoords(spe))[, 1:2, drop = FALSE]
  bin <- dbscan::kNN(xy_node, k = 1, query = xy_cell)$id[, 1]
  gnode <- match(model_data$members, dens_dat$node)

  # Gaussian kernel truncated at four bandwidths in x and in y, applied as
  # a horizontal and a vertical pass
  ra <- floor(4 * bw / lat$sa + 1e-9)
  rr <- floor(4 * bw / lat$sr + 1e-9)
  gx <- exp(-((-ra:ra) * lat$sa)^2 / (2 * bw^2))
  gy <- exp(-((-rr:rr) * lat$sr)^2 / (2 * bw^2))

  # problem for the C++ kernel, with 0 based indices; within each group the
  # cells of the cell types of interest come first, ordered by cell type
  o <- order(gid, is.na(ct_idx), ct_idx)
  ngroup <- max(gid)
  ct_o <- ct_idx[o]
  keep <- !is.na(ct_o)
  overall <- if ("density_overall" %in% colnames(model_data)) {
    as.double(model_data$density_overall)
  }
  prob <- list(
    goff = c(0L, cumsum(tabulate(gid, nbins = ngroup))),
    labs = ct_o[keep] - 1L,
    loff = c(0L, cumsum(tabulate(gid[o][keep], nbins = ngroup))),
    cell_a = lat$a[bin[o]], cell_r = lat$r[bin[o]],
    grid_a = lat$a[gnode], grid_r = lat$r[gnode],
    lattice = as.integer(c(max(lat$a), max(lat$r), ra, rr, lat$cs,
                           lat$parity)),
    gx = gx, gy = gy, overall = overall,
    roi = comp_int - 1L, nroi = nCpnts, nct = nCT)

  pairs <- utils::combn(nCT, 2)
  nP <- ncol(pairs)
  ip <- seq_len(nP)
  ir <- nP + seq_len(nP * nCpnts)
  cpu_threads <- max(1L, as.integer(cpu_threads))

  # observed statistics (the pooled one of every pair, then every pair per
  # ROI) and densities
  res_obs <- .Call("C_corDensityPerm", prob, 0, 0, 0L, 0L, 0, cpu_threads)
  obs_all <- res_obs[[1]]
  D_obs <- res_obs[[2]]
  stored <- as.matrix(model_data[, den_cols, drop = FALSE])
  okr <- stats::complete.cases(stored)
  rc <- suppressWarnings(min(diag(stats::cor(D_obs[okr, , drop = FALSE],
                                             stored[okr, , drop = FALSE]))))
  if (is.finite(rc) && rc < 0.99) {
    warning("Densities rebuilt for the permutation test differ from those ",
            "stored by gridDensity() (minimum correlation ", signif(rc, 3),
            "). The test is internally consistent but does not use the ",
            "stored densities exactly.")
  }
  obs <- list(pooled = obs_all[ip], roi = matrix(obs_all[ir], nP, nCpnts))

  # permutations run in chunks (a multiple of the kernel's block size of 16,
  # so the result does not depend on the chunking or the number of threads)
  seed <- floor(stats::runif(1) * 2^31)
  acc <- matrix(0, length(obs_all), 4L)       # ge, le, n, sum
  chunk <- 1024L
  for (b0 in seq(0L, nperm - 1L, by = chunk)) {
    nb <- min(chunk, nperm - b0)
    acc <- .Call("C_corDensityPerm", prob, obs_all, acc, as.integer(b0),
                 as.integer(nb), seed, cpu_threads)
    if (trace) cat("permutation", b0 + nb, "of", nperm, "\n")
  }
  ge_p <- acc[ip, 1]; le_p <- acc[ip, 2]; n_p <- acc[ip, 3]; sum_p <- acc[ip, 4]
  ge_r <- matrix(acc[ir, 1], nP, nCpnts); le_r <- matrix(acc[ir, 2], nP, nCpnts)
  n_r <- matrix(acc[ir, 3], nP, nCpnts); sum_r <- matrix(acc[ir, 4], nP, nCpnts)

  pval <- function(cnt, n, o) ifelse(is.na(o) | n == 0, NA_real_, (1 + cnt) / (1 + n))
  expct <- function(s, n) ifelse(n == 0, NA_real_, s / n)
  nm <- vapply(den_cols, .ctName, character(1), USE.NAMES = FALSE)

  exp_r <- expct(sum_r, n_r)
  pPos_r <- pval(ge_r, n_r, obs$roi)
  pNeg_r <- pval(le_r, n_r, obs$roi)
  result.ROI <- data.frame(
    celltype1 = rep(nm[pairs[1, ]], each = nCpnts),
    celltype2 = rep(nm[pairs[2, ]], each = nCpnts),
    ROI = rep(cpnts, times = nP),
    ngrid = rep(comp_n, times = nP),
    cor.coef = as.vector(t(obs$roi)),
    cor.expected = as.vector(t(exp_r)),
    shift = as.vector(t(obs$roi - exp_r)),
    p.Pos = as.vector(t(pPos_r)),
    p.Neg = as.vector(t(pNeg_r))
  )

  # consistency of the direction across ROIs, per pair (ROIs with a defined
  # correlation only; significance at the 0.05 level, one sided)
  valid <- !is.na(obs$roi) & !is.na(exp_r)
  cnt <- function(m) as.integer(rowSums(valid & !is.na(m) & m))
  exp_p <- expct(sum_p, n_p)
  result.overall <- data.frame(
    celltype1 = nm[pairs[1, ]],
    celltype2 = nm[pairs[2, ]],
    cor.coef = obs$pooled,
    cor.expected = exp_p,
    shift = obs$pooled - exp_p,
    p.Pos = pval(ge_p, n_p, obs$pooled),
    p.Neg = pval(le_p, n_p, obs$pooled),
    n.roi = as.integer(rowSums(valid)),
    n.above = cnt(obs$roi > exp_r),
    n.below = cnt(obs$roi < exp_r),
    n.sig.Pos = cnt(pPos_r < 0.05),
    n.sig.Neg = cnt(pNeg_r < 0.05)
  )
  list(ROI = result.ROI, overall = result.overall)
}
