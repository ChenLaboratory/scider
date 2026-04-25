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
#' @param pt.shape shape of points.

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
                   pt.shape = 16,
                   pt.size = 1,
                   pt.alpha = 0.6,
                   label = NULL,
                   xlab = NULL,
                   ylab = NULL,
                   cols.scale=NULL) {
  if(!length(SingleCellExperiment::reducedDim(spe))) {
    stop("No dimensionality reduction found.")
  }
  dimred <- dimred %||% SingleCellExperiment::reducedDimNames(spe)[[1]]
  
  toplot <- as.data.frame(SingleCellExperiment::reducedDim(spe,dimred)[,dims])
  colnames(toplot) <- c("x", "y")
  cdata <- as.data.frame(SummarizedExperiment::colData(spe))
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
  
  # Colours
  if (!is.null(group) && is.null(col.p)) {
    n_colour <- length(unique(group))
    if (is.null(cols)) { # Default palette
      col.p <- `if`(isContinuous, col.spec, selectColor(n_colour))
    } else if (is.function(cols)) { # cols is function
      col.p <- cols(n_colour)
    } else { # cols is vector
      col.p <- `if`(isContinuous, cols, rep_len(cols,n_colour))
    }
  }

  #labels
  xlab <- xlab %||% paste(dimred,dims[1])
  ylab <- ylab %||% paste(dimred,dims[2])
  
  # !!group prevents name-clashing in case toplot also has a 'group' column
  p <- ggplot2::ggplot(toplot,aes(x=x, y=y, color=!!group)) +
    ggplot2::geom_point(
      shape = pt.shape,
      size = pt.size,
      alpha = pt.alpha,
    ) +
    labs(x = xlab, y = ylab, color = label) +
    theme_classic()
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