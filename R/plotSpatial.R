#' Plot cells based on spatial coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param reverseY Reverse y coordinates.
#' @param n Integer value. The number of distinct color to be generated,
#' default is 30.
#' @param pt.color color of points. Must be in colData of spe.
#' @param pt.shape shape of points.
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
#' plotSpatial(spe, pt.color = "cell_type", pt.size = 0.3, pt.alpha = 0.2)
#'
plotSpatial <- function(spe, reverseY = FALSE, n = 30, 
                         pt.color = NULL,
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

  col.p <- selectColor(n)

  #This stop "Coordinate system already present..." warning by coord_fixed()
  cf = coord_fixed()
  cf$default = TRUE

  if (!is.null(pt.color)) {
    p = ggplot2::ggplot(toplot,aes(x=x, y=y,color=.data[[pt.color]]))
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
