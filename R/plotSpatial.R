#' Plot cells based on spatial coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param group.by values to group points by. Must be in colData of spe. 
#' If NULL, will try with 'cols' if available.
#' @param feature Feature(s) to colour points by; must be in rownames(spe). If a
#' vector of more than one feature is supplied, one panel is drawn per feature
#' (see \code{per.scale}, \code{ncol}), as in \code{Seurat::FeaturePlot}.
#' @param assay Name of assay to use for plotting feature.
#' @param type Transformation to apply for the group/feature. Options are "raw"
#' , "log", "cpm", "logcpm", or a function that accepts and returns a vector of
#' the same length. For feature plots the legend title defaults to the matching
#' unit ("Counts", "log2 Cts", "CPM", "log2-CPM"); for a single feature the gene
#' name becomes the plot title. Override the legend with \code{label}.
#' @param cols Colour palette. Can be a vector of colours or a function
#' that accepts an integer n and return n colours.
#' @param highlight Optional cells to emphasise, given as either a vector of
#' group.by levels (characters or cluster numbers), or a logical vector of length
#' ncol(spe) selecting cells directly.
#' Highlighted cells are drawn last (on top) at pt.size.highlight; all other cells
#' are light grey at pt.size.
#' @param cols.highlight Colour(s) for the highlighted cells. Defaults to NULL,
#' which keeps each level's usual group.by palette colour. A single colour (e.g. "red")
#' colours all highlighted cells the same; a vector matching the number of 'highlight'
#' entries gives one colour per level (matched by position).
#' @param pt.shape shape of points.
#' @param pt.size.highlight size of highlighted points (see highlight).
#' @param pt.size size of points.
#' @param pt.alpha alpha of points between 0 and 1.
#' @param label label for the legend
#' @param cols.scale vector of position for color if colors should not be
#' evenly positioned. See \link[ggplot2]{scale_color_gradientn}. Only applicable for continuous values.
#' @param transform Name of a transformation for the continuous colour scale
#' (e.g. "log10", "log1p", "pseudo_log", "sqrt"), passed to
#' \link[ggplot2]{scale_color_gradientn}; the colour spectrum is spaced by the
#' transform while the legend stays in original units. Default "identity" (no
#' transformation). Use "log10" for nicely log-spaced legend breaks (needs
#' positive values); "log1p"/"pseudo_log" tolerate zeros but keep linear breaks.
#' Only affects continuous (feature or numeric group.by) colouring. Use "log10" for log-spaced legend breaks (requires positive
#' values); "log1p"/"pseudo_log" tolerate zeros but keep linear breaks.
#' @param reverseY Logical. Whether to reverse Y coordinates. Default is TRUE
#' if the spe contains an image (even if not plotted) and FALSE if otherwise.
#' @param image Logical. Whether to draw the background tissue image (forwarded
#' to \link[scider]{plotImage}). Default FALSE.
#' @param ncol Number of columns when plotting multiple features. Passed to
#' \link[ggplot2]{facet_wrap} (shared scale) or \link[patchwork]{wrap_plots}
#' (per-panel scales). Default NULL lets the layout be chosen automatically.
#' @param per.scale Logical. For multiple features, whether each panel gets its
#' own colour scale (TRUE, default; like \code{Seurat::FeaturePlot}, via the
#' 'patchwork' package; the image is drawn in every panel) or a single shared
#' colour scale across panels (FALSE).
#' @param range Numeric length-2 (lower, upper) cap for continuous colour values
#' - a single feature, or multiple features with a shared scale; values outside
#' are clamped. A single value is taken as the upper bound. Default NULL (no
#' capping). Ignored when \code{per.scale = TRUE} (each panel auto-scales).
#' @param ... Parameters pass to plotImage
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' plotSpatial(spe, group.by = "cell_type", pt.size = 0.5, pt.alpha = 0.6)
#'
plotSpatial <- function(spe,
                        group.by = NULL,
                        feature = NULL,
                        assay = "counts",
                        type = c("log","raw","cpm","logcpm"),
                        cols = NULL,
                        highlight = NULL,
                        cols.highlight = NULL,
                        pt.shape = 16,
                        pt.size = 0.3,
                        pt.size.highlight = 1,
                        pt.alpha = 1,
                        label = NULL,
                        cols.scale = NULL,
                        reverseY = NULL,
                        image = FALSE,
                        ncol = NULL,
                        per.scale = TRUE,
                        range = NULL,
                        transform = "identity",
                        ...) {
  toplot <- SpatialExperiment::spatialCoords(spe)
  colnames(toplot)[1:2] <- c("x", "y")

  cdata <- SummarizedExperiment::colData(spe)
  
  if ("cell_id" %in% colnames(cdata)) {
    cdata <- cdata[, -which(colnames(cdata) == "cell_id")]
  }
  
  # Resolve the type transform (character preset or user function) and the legend
  # unit label; reused by the single- and multi-feature paths.
  type_chr <- if (is.character(type)) match.arg(type) else NULL
  unit <- if (is.null(type_chr)) "Expr" else .typeLabel(type_chr)
  tfun <- if (is.function(type)) type
          else switch(type_chr,
                      raw    = NULL,
                      log    = function(x) {log2(x + 1)},
                      cpm    = function(x) {(x + 0.5) /
                          colSums(as.matrix(spe@assays@data[[assay]])) * 1e6},
                      logcpm = function(x) {log2((x + 0.5) /
                          colSums(as.matrix(spe@assays@data[[assay]])) * 1e6)})

  # Multiple features: one panel per feature (a la Seurat::FeaturePlot).
  if (!is.null(feature) && length(feature) > 1) {
    miss <- !(feature %in% rownames(spe))
    if (any(miss)) {
      message(paste0(paste(feature[miss], collapse = ", "), " not found. Skipping"))
      feature <- feature[!miss]
    }
    if (length(feature) == 0) stop("None of the features are in rownames(spe).")

    # Per-panel scales: a full plotSpatial() per feature (own scale, own image),
    # combined with patchwork. 'range' does not apply (each panel auto-scales).
    if (per.scale) {
      plots <- lapply(feature, function(f)
        plotSpatial(spe, feature = f, assay = assay, type = type, cols = cols,
                    pt.shape = pt.shape, pt.size = pt.size, pt.alpha = pt.alpha,
                    label = label, cols.scale = cols.scale, reverseY = reverseY,
                    image = image, transform = transform, ...))
      return(patchwork::wrap_plots(plots, ncol = ncol))
    }

    # Single shared colour scale: facet, image repeated across panels.
    exprs <- as.matrix(SummarizedExperiment::assay(spe, assay)[feature, , drop = FALSE])
    nc <- nrow(toplot)
    val <- as.vector(t(exprs))
    if (is.function(tfun)) val <- tfun(val)
    long <- data.frame(
      x = rep(toplot[, "x"], times = length(feature)),
      y = rep(toplot[, "y"], times = length(feature)),
      value = val,
      feature = factor(rep(feature, each = nc), levels = feature))
    if (!is.null(range)) long$value <- .capLimits(long$value, range)
    col.p <- .buildColP(long$value, cols, TRUE)
    p <- plotImage(spe, image = image, reverseY = reverseY, ...) +
      ggplot2::geom_point(data = long, aes(x = x, y = y, color = value),
                          shape = pt.shape, size = pt.size, alpha = pt.alpha) +
      ggplot2::facet_wrap(~ feature, ncol = ncol) +
      labs(x = "x", y = "y", color = label %||% unit) +
      theme_classic() +
      .colorScale(col.p, cols.scale, transform)
    p <- update_bound(p, x = toplot[, "x"], y = toplot[, "y"])
    return(p)
  }

  group <- col.p <- main <- NULL

  # Groups. Order is: colData -> assays -> cols
  if (!is.null(group.by) && group.by %in% colnames(cdata)) {
    group <- cdata[[group.by]]
    if (is.null(label)) label <- group.by
  } else if (!is.null(feature) && feature %in% rownames(spe)) {
    group <- SummarizedExperiment::assay(spe,assay)[feature,]
    # Type transform applies to feature expression only, not to a numeric
    # group.by (e.g. QC metrics or IF intensities).
    if (is.function(tfun)) {
      tryCatch({group <- tfun(group)},
               error = function(e){
                 message("Error when applying 'type'. Skipping 'type'.")
               })
    }
    if (is.null(label)) label <- unit
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

  # Plotting
  pts <- if (is.null(hl_size)) {
    ggplot2::geom_point(
      data = as.data.frame(toplot),
      aes(x=x, y=y, color=!!group), # !! prevent name-clashing if toplot$group exists
      shape = pt.shape, size = pt.size, alpha = pt.alpha)
  } else {
    ggplot2::geom_point(
      data = as.data.frame(toplot),
      aes(x=x, y=y, color=!!group, size=!!hl_size),
      shape = pt.shape, alpha = pt.alpha)
  }
  p <- plotImage(spe, image=image, reverseY=reverseY, ...) + pts +
    labs(x = "x", y = "y", color = label, title = main) +
    theme_classic()
  if (!is.null(hl_size)) p <- p + ggplot2::scale_size_identity()
  p <- update_bound(p,
                    x = toplot[,"x"],
                    y = toplot[,"y"])
  
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

utils::globalVariables(c("x", "y"))