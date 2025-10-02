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
#' @param pt.shape shape of points.
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
                        pt.shape = 16, 
                        pt.size = 0.3, 
                        pt.alpha = 0.5,
                        label = NULL,
                        cols.scale = NULL,
                        reverseY = NULL,
                        ...) {
  toplot <- as.data.frame(SpatialExperiment::spatialCoords(spe))
  colnames(toplot) <- c("x", "y")

  cdata <- as.data.frame(SummarizedExperiment::colData(spe))
  
  if ("cell_id" %in% colnames(cdata)) {
    cdata <- cdata[, -which(colnames(cdata) == "cell_id")]
  }
  
  toplot <- cbind(toplot, cdata) |>
    rownames2col("cell_id")
  
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
  
  # Plotting
  p <- plotImage(spe, reverseY=reverseY, ...) +
    ggplot2::geom_point(
      data = toplot,
      aes(x=x, y=y, color=!!group), # !! prevent name-clashing if toplot$group exists
      shape = pt.shape,
      size = pt.size,
      alpha = pt.alpha,
    ) +
    labs(x = "x", y = "y", color = label) +
    theme_classic()
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