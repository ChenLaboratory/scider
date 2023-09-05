#' Plot model statistics using heatmap.
#'
#' @param model.result A data.frame object.
#' @param stats Character value. Choose either coefficient or t. 
#' Coefficient by default.
#' @param roi Specify the ROI to be plotted. Leave it as default if modelling
#' is not done by per ROI or just want to first ROI to be plotted.
#' @param mode Character vector. Choose either withinROI or acrossROI.
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
#' spe <- corDensity(spe, cal_all = TRUE)
#' 
#' plotCorHeatmap(spe)
#' 
plotCorHeatmap <- function(model.result, 
                        stats = c("cor.coef","t","p.Pos","p.Neg"), 
                        roi = NULL,
                        mode = c("withinROI","acrossROI")){
  
  if (!all(c("cor.coef","p.Pos","p.Neg") %in% colnames(model.result))){
    stop("Please run corDensity before using this function.")
  }
  
  fit_dat <- model.result
  
  if (length(stats) != 1){
    stats <- "cor.coef"
  } else if (!(stats %in% c("cor.coef","t","p.Pos","p.Neg"))){
    stop("stats can only allow either cor.coef, t, p.Pos and p.Neg.")
  }
  
  if(length(mode) == 2L){
    mode <- "withinROI"
  }
  
  if(mode == "withinROI"){
    if ("ROI" %in% colnames(fit_dat)){
      if (is.null(roi)){
        message(paste0("Parameter roi is set to NULL, so ROI #", fit_dat$ROI[1], " is plotted."))
        roi <- fit_dat$ROI[1]
        title <- paste0("Statistics (",stats,") of ROI #", roi)
      } else {
        roi <- roi
        title <- paste0("Statistics (",stats,") of ROI #", roi)
      }
      fit_dat_sub <- fit_dat[fit_dat$ROI == roi,]
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
      as.data.frame(optional = TRUE) |> 
      col2rownames("celltype1")
    
    colnames(plot_dat) <- gsub(paste0("^",stats,"\\."), "", colnames(plot_dat))
    
    plot_dat <- plot_dat[colnames(plot_dat),]
    
    filled_data <- fill_symmetric_na(plot_dat)
  } else if (mode == "acrossROI"){
    fit_dat <- as.data.frame(fit_dat[,c("celltype1","celltype2","ROI",stats)])
    
    fit_dat[["pair"]] <- paste0(fit_dat[["celltype1"]],"-",fit_dat[["celltype2"]])
    
    fit_dat <- fit_dat[,c("pair","ROI",stats)]
    
    filled_data <- reshape(fit_dat, idvar = "pair", timevar = "ROI", 
                           direction = "wide") |>
      col2rownames("pair")
    
    colnames(filled_data) <- gsub(paste0("^",stats,"\\."),"ROI#",colnames(x1))
    
  }
  
  
  
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

