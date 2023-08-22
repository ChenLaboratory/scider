#' Plot cells based on spatial coordinates.
#'
#' @param spe A SpatialExperiment object.
#' @param reverseY Reverse y coordinates.
#' @param n Integer value. The number of distinct color to be generated, default is 30.
#' @param ... Aesthetic mappings to pass to `ggplot2::aes_string()`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' 
#' data("xenium_bc_spe")
#'
#' plotSpatial(spe, shape = ".", color = cell_type, size = 0.3, alpha = 0.2)
#' 
plotSpatial <- function(spe, reverseY = FALSE, n = 30, ...){
  
  toplot <- as.data.frame(SpatialExperiment::spatialCoords(spe))
  
  colnames(toplot) <- c("x", "y")
  
  cdata <- as.data.frame(SummarizedExperiment::colData(spe))
  
  if("cell_id" %in% colnames(cdata)){
    cdata <- dplyr::select(cdata, -cell_id)
  }
  
  toplot <- cbind(toplot, cdata) |>
    tibble::rownames_to_column("cell_id")
  
  aesmap <- rlang::enquos(...)
  
  # split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes <- vapply(aesmap, rlang::quo_is_symbolic, FUN.VALUE = logical(1))
    defaultmap <- lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap <- aesmap[is_aes]
  } else {
    defaultmap <- list()
  }
  
  if (reverseY){
    y_tmp <- toplot[,"y"]
    mid_y <- (max(y_tmp) + min(y_tmp))/2
    final_y <- 2 * mid_y - y_tmp
    toplot[,"y"] <- final_y
  }
  
  p <- ggplot2::ggplot(toplot, aes(x = x, y = y, !!!aesmap))
  
  set.seed(100)
  col.p <- randomcoloR::distinctColorPalette(n)
  
  p <- p +
    do.call(ggplot2::geom_point, defaultmap) +
    scale_color_manual(values=col.p) +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(shape=16, size=5)))
  
  return(p)
}


