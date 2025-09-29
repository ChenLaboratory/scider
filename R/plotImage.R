#' Plot background image of spe
#' @param spe A SpatialExperiment object.
#' @param image Logical. Whether to plot image (if present).
#' @param image_id,sample_id sample and image identifiers for scaling factor. 
#' See \link[SpatialExperiment]{scaleFactors} for default behavior.
#' @param reverseY Logical. Whether to reverse Y coordinates. Default is TRUE 
#' if the spe contains an image (even if not plotted) and FALSE if otherwise.
#' @param crop Whether to crop the plot to the spots.
#' @return a ggplot object if there is a valid image, else NULL.
#' @export
plotImage <- function(spe,
                      image = TRUE,
                      image_id = NULL,
                      sample_id = NULL,
                      reverseY = NULL,
                      crop=TRUE) {
  img <- SpatialExperiment::imgRaster(spe,sample_id,image_id)
  reverseY <- reverseY %||% !is.null(img)
  # No image
  if(is.null(img) || !image) {
    if (!(missing(image_id)&&missing(sample_id))) message("image not found")
    p <- ggplot2::ggplot()+ggplot2::coord_fixed()
    if (reverseY) p <- p + ggplot2::scale_y_reverse()
    return(p)
  }
  
  # Scaling image
  scale <- SpatialExperiment::scaleFactors(spe,sample_id,image_id)
  xlim <- c(0,ncol(img)/scale)
  ylim <- c(0,nrow(img)/scale)
  p <- ggplot2::ggplot() 
  if (reverseY) {
    p <- p + 
      ggplot2::annotation_custom(grid::rasterGrob(img),
                                 xmin = xlim[1], xmax = xlim[2],
                                 ymin = -ylim[1], ymax = -ylim[2]) + 
      ggplot2::scale_y_reverse()
  } else {
    p <- p + 
      ggplot2::annotation_custom(grid::rasterGrob(img),
                                 xmin = xlim[1], xmax = xlim[2],
                                 ymin = ylim[1], ymax = ylim[2])
  }
  
  if (crop) {
    xlim <- range(SpatialExperiment::spatialCoords(spe)[,1])
    ylim<-range(SpatialExperiment::spatialCoords(spe)[,2])
    if (reverseY) ylim = rev(ylim)
  }
  p <- p + ggplot2::coord_fixed(xlim=xlim, ylim=ylim)
  return(p)
}

#' Update the x,y limits of a plot
#' @param p A ggplot() object
#' @param x,y Vectors of new x and y limits.
#' @return a ggplot object
update_bound <- function(p,x=NULL,y=NULL) {
  # Whether reverseY has been applied to p
  reverseY <- any(sapply(p$scales$scales, function(s) {
    inherits(s, "ScaleContinuousPosition") &&
      "y" %in% s$aesthetics &&
      "reverse" %in% s$trans$name
  }))
  p$coordinates$limits$x = range(p$coordinates$limits$x,x)
  p$coordinates$limits$y = range(p$coordinates$limits$y,y)
  
  if (reverseY) p$coordinates$limits$y <- rev(p$coordinates$limits$y)
  return(p)
}
