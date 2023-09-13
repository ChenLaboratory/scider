#' Plot contour lines.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of length 1 of the cell type of interest (COIs).
#' @param overlay Character vector. Either plot overlay on density or cell. 
#' By default is cell.
#' @param id A character. The name of the column of colData(spe) containing the cell type identifiers.
#' Set to cell_type by default.
#' @param sub.level Character vector. Subset on specific level.
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
                        overlay = c("cell","density"),
                        id = "cell_type", 
                        sub.level = NULL, ...){
  
  if (length(coi) > 1L) {
    stop("coi must be of length 1!")
  }
  
  if (!(coi %in% colData(spe)[[id]])){
    stop("coi not in colData(spe)[[id]]!")
  }
  
  coi_clean <- janitor::make_clean_names(coi)
  coi_clean_contour <- paste(coi_clean, "contour", sep = "_")
  
  if(!coi_clean_contour %in% names(spe@metadata)){
    stop("Contour of coi doesn't exist. Please run getContour() first!")
  }

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
    p <- plotDensity(spe, coi = coi, ...)
  } else {
    stop("Overlay should either be cell or density.")
  }
  
  col.p <- grDevices::colorRampPalette(col.spec)(length(unique(contour_data$level)))
  
  if(is.null(sub.level)){
    suppressMessages(
      p <- p +
        ggplot2::geom_path(data = contour_data, 
                           ggplot2::aes(x = x, y = y, group = group, 
                                        color = level)) +
        scale_color_manual(name = "Density level", values = rev(col.p))
    )
  } else {
    if(length(sub.level) == 1L & sub.level %in% contour_data$level){
      suppressMessages(
        p <- p +
          ggplot2::geom_path(data = contour_data, 
                             ggplot2::aes(x = x, y = y, group = group, 
                                          color = level == sub.level)) +
          scale_color_manual(name = paste0("level",sub.level," density"),
                             values = c("royalblue","tomato2"))
      )
    } else {
      stop("The length sub.level is expected to be 1 and 
           should be included in contour_data$level.")
    }
  }
  
  p <- p +
    theme_classic() +
    labs(x = "x", y = "y") +
    ggtitle(coi) 
  return(p)
}

