#' Plot model statistics using heatmap.
#'
#' @param spe A SpatialExperiment object.
#' @param stats Character value. Choose either coefficient or t. 
#' Coefficient by default.
#' @param roi Specify the ROI to be plotted. Leave it as default if modelling
#' is not done by per ROI or just want to first ROI to be plotted.
#'
#' @return A pheatmap object.
#' @export
#'
#' @examples
#' 
#' data("xenium_bc_spe")
#'
#' coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
#' 
#' spe <- gridDensity(spe, coi = coi)
#' 
#' spe <- findROI(spe, coi = coi, method = "walktrap")
#' 
#' spe <- corrDensity(spe, cal_all = TRUE)
#' 
#' plotModStat(spe)
#' 
plotModStat <- function(spe, stats = c("coefficient","t"), roi = NULL){
  
  if (!("model_result" %in% names(spe@metadata))){
    stop("Please run corrDensity before using this function.")
  }
  
  fit_dat <- as.data.frame(spe@metadata$model_result)
  
  if (length(stats) != 1){
    stats <- "coefficient"
  } else if (!(stats %in% c("coefficient","t"))){
    stop("stats can only allow either coefficient or t.")
  }
  
  if (!(stats %in% colnames(fit_dat))){
    stop(paste0("stats: ", stats, " is not in the colnames of spe@metadata$model_result."))
  }
  
  if ("roi" %in% colnames(fit_dat)){
    if (is.null(roi)){
      message(paste0("Parameter roi is set to NULL, so ROI #", fit_dat$roi[1], " is plotted."))
      roi <- fit_dat$roi[1]
      title <- paste0("Statistics (",stats,") of ROI #", roi)
    } else {
      roi <- roi
      title <- paste0("Statistics (",stats,") of ROI #", roi)
    }
    fit_dat_sub <- fit_dat[fit_dat$roi == roi,]
  } else {
    fit_dat_sub <- fit_dat
    title <- paste0("Statistics (",stats,")")
  }
  
  
  
  plot_dat <- fit_dat_sub[, c("celltype1","celltype2",stats)]
  
  cts <- unique(c(plot_dat$celltype1, plot_dat$celltype2))
  
  sup_dat <- data.frame("celltype1" = cts,
                        "celltype2" = cts)
  
  sup_dat[,stats] <- 0
  
  plot_dat <- rbind(plot_dat, sup_dat)
  
  plot_dat <- reshape(plot_dat, 
                      idvar = "celltype1", 
                      timevar = "celltype2", 
                      direction = "wide") |>
    col2rownames("celltype1")
  
  colnames(plot_dat) <- gsub(paste0("^",stats,"\\."), "", colnames(plot_dat))
  
  plot_dat <- plot_dat[colnames(plot_dat),]
  
  filled_data <- fill_symmetric_na(plot_dat)
  
  paletteLength <- 100
  hmColor <- rev(colorRampPalette(c("#ED254EFF","#EF6079FF", "#F1F4FFFF", "#97B3D0FF", "#011936FF"))(paletteLength))
  #myBreaks <- c(seq(min(filled_data), 0, length.out=ceiling(paletteLength/2) + 1), 
  #              seq(max(filled_data)/paletteLength, max(filled_data), length.out=floor(paletteLength/2)))
  
  pheatmap::pheatmap(filled_data, angle_col = 45, border_color = "white",
                     color = hmColor, #breaks = myBreaks,
                     main = title)
}


fill_symmetric_na <- function(mat) {
  is_na <- is.na(mat)
  mat[is_na] <- t(mat)[is_na]
  return(mat)
}

