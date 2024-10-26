#' Plot cells based on spatial coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param reverseY Reverse y coordinates.
#' @param group.by values to group points by. Must be in colData of spe. 
#' If NULL, will try with 'cols' if available.
#' @param pt.shape shape of points.
#' @param cols Colour palette. Can be a vector of colours or a function 
#' that accepts an integer n and return n colours.
#' @param pt.size size of points.
#' @param pt.alpha alpha of points between 0 and 1.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' plotSpatial(spe, group.by = "cell_type", pt.size = 0.5, pt.alpha = 0.6)
#'
plotSpatial <- function(spe, reverseY = FALSE,
                         group.by = NULL,
                         cols = NULL,
                         pt.shape = 16, 
                         pt.size = 0.3, 
                         pt.alpha = 0.5) {
  toplot <- as.data.frame(SpatialExperiment::spatialCoords(spe))

  colnames(toplot) <- c("x", "y")

  cdata <- as.data.frame(SummarizedExperiment::colData(spe))

  if ("cell_id" %in% colnames(cdata)) {
    cdata <- cdata[, -which(colnames(cdata) == "cell_id")]
  }

  toplot <- cbind(toplot, cdata) |>
    rownames2col("cell_id")

  if (reverseY) {
    y_tmp <- toplot[, "y"]
    mid_y <- (max(y_tmp) + min(y_tmp)) / 2
    final_y <- 2 * mid_y - y_tmp
    toplot[, "y"] <- final_y
  }

  # Groups
  if (!is.null(group.by)) {
    group = toplot[[group.by]]
  } else if (!is.null(cols)) {
    group = factor(rep_len(cols,nrow(toplot)),levels=unique(cols))
  } else {
    group = NULL
  }
  
  isContinuous = !is.null(group.by) && is.numeric(toplot[[group.by]])
  n_colour = length(unique(group))
  if (is.null(cols)&&is.null(group.by)) {
    col.p = NULL
  } else if (is.null(cols)) {
    if (isContinuous) col.p <- col.spec
    else col.p <- selectColor(n_colour)
  } else if (is.function(cols)) {
    col.p <- as.character(cols(n_colour))
  } else {
    col.p <- as.character(cols)
    if (!is.null(group.by) && !isContinuous) col.p <- rep_len(col.p,n_colour)
    else if (!isContinuous) col.p <- rep_len(unique(col.p), n_colour)
  }
  
  #This stop "Coordinate system already present..." warning by coord_fixed()
  cf = coord_fixed()
  cf$default = TRUE

  p = ggplot2::ggplot(toplot,aes(x=x, y=y, color=group)) +
    ggplot2::geom_point(
      shape = pt.shape,
      size = pt.size,
      alpha = pt.alpha,
      ) +
    labs(x = "x", y = "y", color = group.by) +
    theme_classic() +
    cf
  if (isContinuous) {
    p = p + scale_color_gradientn(colours = rev(col.p))
  } else {
    p = p + scale_color_manual(values = col.p) +
      guides(colour = guide_legend(override.aes = list(
        shape = 16,
        size = 5
      )))
  }
  
  return(p)
}

utils::globalVariables(c("x", "y"))
