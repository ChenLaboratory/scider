#' Plot contour lines.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of length 1 of the cell type of interest (COIs).
#' @param overlay Character vector. Either plot overlay on density or cell. 
#' By default is cell.
#' @param sub_level Character vector. Subset on specific level.
#' @param seed Integer. Used for random color selection.
#' @param ... Aesthetic mappings to pass to `ggplot2::aes_string()`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#' 
#' spe <- gridDensity(spe)
#' 
#' coi <- "Breast cancer"
#' 
#' spe <- getContour(spe, coi = coi)
#' 
#' plotContour(spe, coi = coi, size = 0.3, alpha = 0.2)
#' 
plotContour <- function(spe, 
                        coi, 
                        id = "cell_type", 
                        overlay = c("cell","density"),
                        sub_level = NULL, ...){
  
  if (length(coi) > 1L) {
    stop("coi must be of length 1!")
  }
  
  coi_clean <- janitor::make_clean_names(coi)
  coi_clean_contour <- paste(coi_clean, "contour", sep = "_")
  
  contour_data <- as.data.frame(spe@metadata[[coi_clean_contour]])
  levs <- unique(contour_data$level)
  nlevs <- length(levs)
  
  if(length(overlay) == 2){
    overlay <- "cell"
  }
  
  if(overlay == "cell"){
    sub <- grep(coi, colData(spe)[[id]])
    p <- plotSpatial(spe[, sub], ...)
  } else if (overlay == "density"){
    p <- plotDensity(spe, coi = coi)
  } else {
    stop("Overlay should either be cell or density.")
  }
  
  if(is.null(sub_level)){
    suppressMessages(
      p <- p +
        ggplot2::geom_path(data = contour_data, 
                           ggplot2::aes(x = x, y = y, group = group, 
                                        color = level)) +
        labs(color = "Density level")
    )
  } else {
    if(length(sub_level) == 1L & sub_level %in% contour_data$level){
      suppressMessages(
        p <- p +
          ggplot2::geom_path(data = contour_data, 
                             ggplot2::aes(x = x, y = y, group = group, 
                                          color = level == sub_level)) +
          scale_color_manual(name = paste0("level",sub_level," density"),
                             values = c("royalblue","tomato2"))
      )
    } else {
      stop("The length sub_level is expected to be 1 and 
           should be included in contour_data$level.")
    }
  }
  
  p <- p +
    theme_classic() +
    labs(x = "x", y = "y") +
    ggtitle(coi) 
  return(p)
}

