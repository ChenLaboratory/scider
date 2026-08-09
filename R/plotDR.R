#' Plot reduced dimensions.
#'
#' plotDR is the main function for plotting reduced dimension. Others are
#' wrapper functions for convenience.
#' @param spe A SpatialExperiment object.
#' @param dimred Name of the reduced dimension in \link[SingleCellExperiment]{reducedDims} 
#' @param dims Numeric vector length 2 for the dimensions to be plotted. Default to first two dimensions
#' @param group.by values to group points by. Must be in colData of spe. 
#' If NULL, will try with 'cols' if available.
#' @param cols Colour palette. Can be a vector of colours or a function
#' that accepts an integer n and return n colours.
#' @param feature Feature(s) to colour points by; must be in rownames(spe). If a
#' vector of more than one feature is supplied, one reduced-dimension panel is
#' drawn per feature, faceted (see \code{ncol}), with a shared colour scale.
#' @param assay Name of assay to use for plotting feature. Default "counts".
#' @param type Transformation applied to feature expression: "log" (default,
#' log2(1+x)), "raw" (no transform), "cpm", or "logcpm" (log2 CPM). For multiple
#' features the legend title defaults to the matching unit ("log2 Cts", "Counts",
#' "CPM", "log2-CPM"); for a single feature the unit is the legend title and the
#' gene name becomes the plot title. Override the legend with \code{label}.
#' @param ncol Number of columns when plotting multiple features. Passed to
#' \link[ggplot2]{facet_wrap} (shared scale) or \link[patchwork]{wrap_plots}
#' (per-panel scales). Default NULL lets the layout be chosen automatically.
#' @param per.scale Logical. For multiple features, whether each panel gets its
#' own colour scale (TRUE, default; like \code{Seurat::FeaturePlot}, via the
#' 'patchwork' package) or a single shared colour scale across panels (FALSE).
#' @param range Numeric length-2 (lower, upper) cap for continuous colour values
#' - a single feature, or multiple features with a shared scale; values outside
#' are clamped. A single value is taken as the upper bound. Default NULL (no
#' capping). Ignored when \code{per.scale = TRUE} (each panel auto-scales).
#' @param highlight Optional cells to emphasise, given as either a vector of
#' group.by levels (characters or cluster numbers), or a logical vector of length
#' ncol(spe) selecting cells directly. Highlighted cells are drawn last (on top) 
#' at pt.size.highlight; all other cells are light grey at pt.size.
#' @param cols.highlight Colour(s) for the highlighted cells. Defaults to NULL,
#' which keeps each level's usual group.by palette colour. A single colour (e.g. "red")
#' colours all highlighted cells the same; a vector matching the number of 'highlight' 
#' entries gives one colour per level (matched by position).
#' @param pt.shape shape of points.
#' @param pt.size.highlight size of highlighted points (see highlight).
#' @param pt.size size of points.
#' @param pt.alpha alpha of points between 0 and 1.
#' @param label label for the legend
#' @param xlab label for the x-axis
#' @param ylab label for the y-axis
#' @param cols.scale vector of position for color if colors should not be
#' evenly positioned. See \link[ggplot2]{scale_color_gradientn}. Only applicable for continuous values.
#' @param transform Name of a transformation for the continuous colour scale
#' (e.g. "log10", "log1p", "pseudo_log", "sqrt"), passed to
#' \link[ggplot2]{scale_color_gradientn}; the colour spectrum is spaced by the
#' transform while the legend stays in original units. Default "identity" (no
#' transformation). Use "log10" for nicely log-spaced legend breaks (needs
#' positive values); "log1p"/"pseudo_log" tolerate zeros but keep linear breaks.
#' Only affects continuous (feature or numeric group.by) colouring.
#' @param ... Additional arguments pass to plotDR
#' @return A ggplot object.
#' 
#' @rdname plotDR
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#' spe = runUMAP(spe)
#' plotDR(spe, group.by = "cell_type")
#'
plotDR <- function(spe, dimred = NULL,
                   dims = c(1,2),
                   group.by = NULL,
                   feature = NULL,
                   assay = "counts",
                   type = c("log", "raw", "cpm", "logcpm"),
                   cols = NULL,
                   highlight = NULL,
                   cols.highlight = NULL,
                   pt.shape = 16,
                   pt.size = 0.7,
                   pt.size.highlight = 1,
                   pt.alpha = 0.6,
                   label = NULL,
                   xlab = NULL,
                   ylab = NULL,
                   cols.scale=NULL,
                   ncol = NULL,
                   per.scale = TRUE,
                   range = NULL,
                   transform = "identity") {
  if(!length(rds <- SingleCellExperiment::reducedDimNames(spe))) {
    stop("No dimensionality reduction found.")
  }
  dimred <- dimred %||% rds[[1]]
  type <- match.arg(type)

  toplot <- SingleCellExperiment::reducedDim(spe,dimred)[,dims]
  colnames(toplot) <- c("x", "y")
  # Use a data.frame so `toplot$x`/`toplot$y` work (reducedDim returns a matrix,
  # on which `$` errors); matrix-style indexing below still works on a df.
  toplot <- as.data.frame(toplot)
  cdata <- SummarizedExperiment::colData(spe)

  # Multiple features: one panel per feature (a la Seurat::FeaturePlot).
  if (!is.null(feature) && length(feature) > 1) {
    miss <- !(feature %in% rownames(spe))
    if (any(miss)) {
      message(paste0(paste(feature[miss], collapse = ", "), " not found. Skipping"))
      feature <- feature[!miss]
    }
    if (length(feature) == 0) stop("None of the features are in rownames(spe).")
    exprs <- as.matrix(SummarizedExperiment::assay(spe, assay)[feature, , drop = FALSE])
    libsize <- .libSize(spe, assay, type)
    unit <- label %||% .typeLabel(type)
    xlab <- xlab %||% paste(dimred, dims[1])
    ylab <- ylab %||% paste(dimred, dims[2])

    # Per-panel scales: build a separate plot (own colour scale) per feature and
    # combine with patchwork. 'range' does not apply (each panel auto-scales).
    if (per.scale) {
      plots <- lapply(feature, function(f) {
        v <- .transformExpr(exprs[f, ], type, libsize)
        d <- data.frame(x = toplot$x, y = toplot$y, value = v)
        ggplot2::ggplot(d, aes(x = x, y = y, color = value)) +
          ggplot2::geom_point(shape = pt.shape, size = pt.size, alpha = pt.alpha) +
          labs(x = xlab, y = ylab, color = unit, title = f) +
          theme_classic() +
          .colorScale(.buildColP(v, cols, TRUE), cols.scale, transform)
      })
      return(patchwork::wrap_plots(plots, ncol = ncol))
    }

    # Single shared colour scale: facet. 'range' optionally caps the values.
    nc <- nrow(toplot)
    ls.long <- if (is.null(libsize)) NULL else rep(libsize, times = length(feature))
    long <- data.frame(
      x = rep(toplot$x, times = length(feature)),
      y = rep(toplot$y, times = length(feature)),
      value = .transformExpr(as.vector(t(exprs)), type, ls.long),
      feature = factor(rep(feature, each = nc), levels = feature))
    if (!is.null(range)) long$value <- .capLimits(long$value, range)
    p <- ggplot2::ggplot(long, aes(x = x, y = y, color = value)) +
      ggplot2::geom_point(shape = pt.shape, size = pt.size, alpha = pt.alpha) +
      ggplot2::facet_wrap(~ feature, ncol = ncol) +
      labs(x = xlab, y = ylab, color = unit) +
      theme_classic() +
      .colorScale(.buildColP(long$value, cols, TRUE), cols.scale, transform)
    return(p)
  }

  group <- col.p <- main <- NULL
  # Groups. Order is: colData -> assays -> cols
  if (!is.null(group.by) && group.by %in% colnames(cdata)) {
    group <- cdata[[group.by]]
    if (is.null(label)) label <- group.by
  } else if (!is.null(feature) && feature %in% rownames(spe)) {
    group <- SummarizedExperiment::assay(spe,assay)[feature,]
    group <- .transformExpr(group, type, .libSize(spe, assay, type))
    if (is.null(label)) label <- .typeLabel(type)
    main <- feature
  } else if (!is.null(cols) && !is.function(cols)) {
    group <- factor(rep_len(cols,nrow(toplot)),levels=unique(cols))
    col.p <- rep_len(unique(cols), length(unique(cols)))
  }
  isContinuous <- is.numeric(group)
  if (isContinuous && !is.null(range)) group <- .capLimits(group, range)

  # Highlight a subset of groups, or build the normal palette.
  # ('unassigned' cells are coloured black by default)
  hl_size <- NULL
  if (!is.null(highlight) && !isContinuous) {
    hl <- .prepHighlight(group, highlight, cols, cols.highlight,
                         pt.size, pt.size.highlight)
    group   <- hl$group[hl$order]
    col.p   <- hl$col.p
    hl_size <- hl$size[hl$order]
    toplot  <- toplot[hl$order, , drop = FALSE]
  } else if (!is.null(group) && is.null(col.p)) {
    col.p <- .buildColP(group, cols, isContinuous)
  }

  #labels
  xlab <- xlab %||% paste(dimred,dims[1])
  ylab <- ylab %||% paste(dimred,dims[2])

  # !!group prevents name-clashing in case toplot also has a 'group' column
  p <- ggplot2::ggplot(as.data.frame(toplot),
                       aes(x=x, y=y, color=!!group))
  if (is.null(hl_size)) {
    p <- p + ggplot2::geom_point(shape = pt.shape, size = pt.size, alpha = pt.alpha)
  } else {
    p <- p + ggplot2::geom_point(aes(size = !!hl_size), shape = pt.shape,
                                 alpha = pt.alpha) +
      ggplot2::scale_size_identity()
  }
  p <- p + labs(x = xlab, y = ylab, color = label, title = main) + theme_classic()
  if (isContinuous) {
    p <- p + .colorScale(col.p, cols.scale, transform)
  } else {
    p <- p + scale_color_manual(values = col.p) +
      guides(colour = guide_legend(override.aes = list(
        shape = 16,
        size = 5
      )))
  }
  return(p)
}

#' @rdname plotDR
#' @aliases plotUMAP
#' 
#' @export
plotUMAP <- function (spe,dimred="UMAP",...) {
  plotDR(spe,dimred,...)
}

#' @rdname plotDR
#' @aliases plotPCA
#' 
#' @export
plotPCA <- function (spe,dimred="PCA",...) {
  args <- list(...)
  args$dims = args$dims %||% c(1,2)
  percentVar = attr(SingleCellExperiment::reducedDim(spe,dimred),"percentVar")[args$dims]
  if (!is.null(percentVar)) {
    percentVar = round(percentVar)
    args$xlab <- args$xlab %||% paste0(dimred," ",args$dims[1]," (",percentVar[1],"%)")
    args$ylab <- args$ylab %||% paste0(dimred," ",args$dims[2]," (",percentVar[2],"%)")
  }
  
  args$spe = spe
  args$dimred = dimred
  do.call(plotDR,args)
}


# Library sizes (per cell) for cpm/logcpm transforms; NULL when not needed.
.libSize <- function(spe, assay, type) {
  if (type %in% c("cpm", "logcpm"))
    colSums(as.matrix(spe@assays@data[[assay]]))
  else NULL
}

# Transform an expression vector by 'type'. libsize must align with x for
# cpm/logcpm (recycled per cell); may be NULL otherwise. "log" is log2(1+x),
# matching plotSpatial().
.transformExpr <- function(x, type, libsize = NULL) {
  switch(type,
         raw    = x,
         log    = log2(x + 1),
         cpm    = (x + 0.5) / libsize * 1e6,
         logcpm = log2((x + 0.5) / libsize * 1e6))
}

# Continuous colour scale, handling the ggplot2 3.5 rename of the transform
# argument (trans -> transform).
.colorScale <- function(col.p, values = NULL, trans = "identity") {
  args <- list(colours = rev(col.p), values = values)
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    args$transform <- trans
  else
    args$trans <- trans
  do.call(ggplot2::scale_color_gradientn, args)
}

# Legend unit label for each 'type'.
.typeLabel <- function(type) {
  switch(type,
         raw    = "Counts",
         log    = "log2 Cts",
         cpm    = "CPM",
         logcpm = "log2-CPM")
}

utils::globalVariables(c("value", "feature"))