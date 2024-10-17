#' Plot cells based on spatial coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param reverseY Reverse y coordinates.
#' @param n Integer value. The number of distinct colour to be generated,
#' default is 30.
#' @param colour.by values to colour points by. Must be in colData of spe.
#' @param pt.shape shape of points.
#' @param pt.colour Colour palette. Can be a vector of colours or a function 
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
#' plotSpatial(spe, colour.by = "cell_type", pt.size = 0.3, pt.alpha = 0.2)
#'
plotSpatial <- function(spe, reverseY = FALSE, n = 30, 
                         colour.by = NULL,
                         pt.colour = NULL,
                         pt.shape = 16, 
                         pt.size = 0.3, 
                         pt.alpha = 0.2) {
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

  col.p = NULL
  if (!is.null(colour.by)) {
    n_colour = length(unique(toplot[[colour.by]]))
    if (is.null(pt.colour)) {
      col.p <- selectColor(n_colour)
    } else if (is.function(pt.colour)) {
      col.p <- pt.colour(n_colour)
    } else {
      col.p <- rep_len(pt.colour, n_colour)
    }
  }

  #This stop "Coordinate system already present..." warning by coord_fixed()
  cf = coord_fixed()
  cf$default = TRUE

  if (!is.null(colour.by)) {
    p = ggplot2::ggplot(toplot,aes(x=x, y=y,color=.data[[colour.by]]))
  } else {
    p = ggplot2::ggplot(toplot,aes(x=x, y=y))
  }
  p = p +
    ggplot2::geom_point(
      shape = pt.shape,
      size = pt.size,
      alpha = pt.alpha
      ) +
    scale_color_manual(values = col.p) +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(
      shape = 16,
      size = 5
    ))) +
    cf

  return(p)
}

utils::globalVariables(c("x", "y"))
