#' Plot cells based on spatial coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param group.by values to group points by. Must be in colData of spe. 
#' If NULL, will try with 'cols' if available.
#' @param feature Feature to group polygons by. Must be in rownames(spe).
#' @param assay Name of assay to use for plotting feature.
#' @param type Transformation to apply for the group/feature. Options are "raw"
#' , "log", "cpm", "logcpm", or a function that accepts and returns a vector of 
#' the same length.
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
#' @param reverseY Logical. Whether to reverse Y coordinates. Default is TRUE 
#' if the spe contains an image (even if not plotted) and FALSE if otherwise.
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
                        type = c("raw","log","cpm","logcpm"),
                        cols = NULL,
                        highlight = NULL,
                        cols.highlight = NULL,
                        pt.shape = 16,
                        pt.size = 0.3,
                        pt.size.highlight = 1,
                        pt.alpha = 0.5,
                        label = NULL,
                        cols.scale = NULL,
                        reverseY = NULL,
                        ...) {
  toplot <- SpatialExperiment::spatialCoords(spe)
  colnames(toplot)[1:2] <- c("x", "y")

  cdata <- SummarizedExperiment::colData(spe)
  
  if ("cell_id" %in% colnames(cdata)) {
    cdata <- cdata[, -which(colnames(cdata) == "cell_id")]
  }
  
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
  
  # Type
  if (is.character(type)) {
    type <- switch(match.arg(type),
                   raw = NULL,
                   log = function(x) {log2(x+1)},
                   cpm = function(x) {
                     (x+0.5)/colSums(as.matrix(spe@assays@data[[assay]]))*1e6
                   },
                   logcpm = function(x) {
                     log2((x+0.5)/colSums(as.matrix(spe@assays@data[[assay]]))*1e6)
                   })
  }
  if(is.function(type)) {
    tryCatch({group <- type(group)},
             error = function(e){
               message("Error when applying 'type'. Skipping 'type'.")
             })
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
  p <- plotImage(spe, reverseY=reverseY, ...) + pts +
    labs(x = "x", y = "y", color = label) +
    theme_classic()
  if (!is.null(hl_size)) p <- p + ggplot2::scale_size_identity()
  p <- update_bound(p,
                    x = toplot[,"x"],
                    y = toplot[,"y"])
  
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

utils::globalVariables(c("x", "y"))