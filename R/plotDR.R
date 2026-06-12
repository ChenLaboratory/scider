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
#' @param feature Feature to group points by. Must be in rownames(spe).
#' @param assay Name of assay to use for plotting feature.
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
                   cols.scale=NULL) {
  if(!length(rds <- SingleCellExperiment::reducedDimNames(spe))) {
    stop("No dimensionality reduction found.")
  }
  dimred <- dimred %||% rds[[1]]
  
  toplot <- SingleCellExperiment::reducedDim(spe,dimred)[,dims]
  colnames(toplot) <- c("x", "y")
  cdata <- SummarizedExperiment::colData(spe)
  toplot <- cbind(toplot, cdata)
  
  group <- col.p <- NULL
  # Groups. Order is: colData -> assays -> cols
  if (!is.null(group.by) && group.by %in% colnames(toplot)) {
    group <- toplot[[group.by]]
    if (is.null(label)) label <- group.by
  } else if (!is.null(feature) && feature %in% rownames(spe)) {
    group <- SummarizedExperiment::assay(spe,assay)[feature,]
    if (is.null(label)) label <- feature
  } else if (!is.null(cols) && !is.function(cols)) {
    group <- factor(rep_len(cols,nrow(toplot)),levels=unique(cols))
    col.p <- rep_len(unique(cols), length(unique(cols)))
  }
  isContinuous <- is.numeric(group)

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
  p <- p + labs(x = xlab, y = ylab, color = label) + theme_classic()
  if (isContinuous) {
    p <- p + scale_color_gradientn(colours = rev(col.p), values = cols.scale)
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