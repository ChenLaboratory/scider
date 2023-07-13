#' Plot cells based on cell position on tissue.
#'
#' @param spe SpatialExperiment object.
#' @param targetcell Optional. Can input ONE specific cell id to zoom-in on the region of a specific cell.
#' @param k_near Optional. If targetcell is specified, the k_near cells around the targetcell will be plotted.
#' @param targetsize Dot size of the targetcell.
#' @param targetshape Shape of the targetcell.
#' @param targetcolor Colour of the targetcell.
#' @param scaleFactor Scale factor to align with the image.
#' @param reverseY Reverse y coordinates.
#' @param plotImg Logical. Set to TRUE to plot the tissue image at the bottom layer.
#' @param sampleid Sample id link to the tissue image.
#' @param image_id Image id link to the tissue image.
#' @param ... aesthetic mappings to pass to `ggplot2::aes_string()`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' plot_tissue(spe, shape = ".")
#'
plot_tissue <- function(spe, targetcell = FALSE, k_near = 100, targetsize = 3,
                        targetshape = 1, targetcolor = "red",
                        scaleFactor = 1, reverseY = TRUE, plotImg = FALSE,
                        sampleid = "V1", image_id = "V1_IDC", ...){
  if (!is(targetcell, "logical")){
    if (length(targetcell) == 1){
      ncells <- find_near_cells(spe, k = k_near, targetCell = targetcell, reportCellID = TRUE,
                                reportDist = FALSE) |>
        unlist()
      
      spe <- spe[,c(targetcell,ncells)]
    }
  }
  
  toplot <- SpatialExperiment::spatialCoords(spe) |>
    as.data.frame()
  
  colnames(toplot) <- c("x", "y")
  
  cdata <- as.data.frame(SummarizedExperiment::colData(spe))
  
  if("cell_id" %in% colnames(cdata)){
    cdata <- dplyr::select(cdata, -cell_id)
  }
  
  toplot <- toplot |>
    cbind(cdata) |>
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
  
  toplot[,"x"] <- scaleFactor * toplot[,"x"]
  toplot[,"y"] <- scaleFactor * toplot[,"y"]
  
  if (reverseY){
    y_tmp <- toplot[,"y"]
    mid_y <- (max(y_tmp) + min(y_tmp))/2
    final_y <- 2 * mid_y - y_tmp
    toplot[,"y"] <- final_y
  }
  
  p <- ggplot2::ggplot(toplot, aes(x = x, y = y, !!!aesmap))
  
  if (plotImg) {
    img_df <- imgData(spe)
    rownames(img_df) <- img_df$sample_id
    
    spi <- img_df[sampleid, "data"]
    img <- imgRaster(spi[[1]])
    
    img_layer <- layer(data = data.frame(sample_id = sampleid),
                       inherit.aes = FALSE,
                       stat = "identity",
                       position = "identity",
                       geom = ggplot2::GeomCustomAnn,
                       params = list(
                         grob = grid::rasterGrob(img),
                         xmin = 0, xmax = ncol(img),
                         ymin = 0, ymax = nrow(img))
    )
    
    p <- p + img_layer
    
  }
  
  p <- p +
    do.call(ggplot2::geom_point, defaultmap) +
    tissue_theme()
  
  if(!is(targetcell, "logical")){
    target_df <- toplot |>
      dplyr::filter(cell_id %in% targetcell)
    
    p <- p +
      ggplot2::geom_point(data = target_df, aes(x, y), size = targetsize, shape = targetshape, color = targetcolor)
  }
  
  return(p)
}



tissue_theme <- function(textScale = 1.1){
  stopifnot(textScale > 0)
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      legend.text = element_text(size = rel(textScale)),
      legend.title = element_text(size = rel(textScale), face = "italic"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}


utils::globalVariables(c("x", "y","cell_id"))

