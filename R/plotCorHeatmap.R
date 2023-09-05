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
#' model_result <- corDensity(spe, cal_all = TRUE)
#' 
#' plotCorHeatmap(model_result)
#' 
plotCorHeatmap <- function(model.result, 
                        stats = c("cor.coef","t","p.Pos","p.Neg"), 
                        roi = "all",
                        cell.type = "all"){
  
  if (!all(c("cor.coef","p.Pos","p.Neg") %in% colnames(model.result))){
    stop("Please run corDensity before using this function.")
  }
  
  fit_dat <- model.result
  
  if (length(stats) != 1){
    stats <- "cor.coef"
  } else if (!(stats %in% c("cor.coef","t","p.Pos","p.Neg"))){
    stop("stats can only allow either cor.coef, t, p.Pos and p.Neg.")
  }
  
  if (all(cell.type != "all")){
    if (length(cell.type) < 2L){
      stop("cell.type must be either all or length larger than 1.")
    }
    if (!all(cell.type %in% unique(c(fit_dat$celltype1, fit_dat$celltype2)))){
      stop("cell.type must be all in celltype1 or celltype2 columns.")
    }
    
    fit_dat <- fit_dat[(fit_dat$celltype1 %in% cell.type) & 
                         (fit_dat$celltype2 %in% cell.type),]
  }
  
  if (!("ROI" %in% colnames(fit_dat))){
    fit_dat <- fit_dat
    title <- paste0("Statistics (",stats,")")
    filled_data <- plot_type1_transformation(fit_dat, stats)
  } else {
    if (all(roi != "all")){
      if (length(roi) == 1L){
        if (!(roi %in% fit_dat[["ROI"]])){
          stop("roi is not in the ROI column.")
        }
        fit_dat <- fit_dat[fit_dat$ROI == roi,]
        title <- paste0("Statistics (",stats,") of ROI #", roi)
        filled_data <- plot_type1_transformation(fit_dat, stats)
      } else {
        if (!all(roi %in% fit_dat[["ROI"]])){
          stop("rois are not all in the ROI column.")
        }
        fit_dat <- fit_dat[fit_dat$ROI %in% roi,]
        title <- paste0("Statistics (",stats,") of ROI #", paste(roi,
                                                                 collapse=","))
        filled_data <- plot_type2_transformation(fit_dat, stats)
      }
    } else {
      title <- paste0("Statistics (",stats,")" )
      filled_data <- plot_type2_transformation(fit_dat, stats)
    }
  }
  
  paletteLength <- 100
  hmColor <- rev(colorRampPalette(c("#ED254EFF","#EF6079FF", 
                                    "#F1F4FFFF", "#97B3D0FF", 
                                    "#011936FF"))(paletteLength))
  #myBreaks <- c(seq(min(filled_data), 0, length.out=ceiling(paletteLength/2)+1), 
  #              seq(max(filled_data)/paletteLength, max(filled_data), 
  #                  length.out=floor(paletteLength/2)))
  
  if(nrow(filled_data) == 1L){
    pheatmap::pheatmap(filled_data, angle_col = 45, border_color = "white",
                       color = hmColor, #breaks = myBreaks,
                       main = title, cluster_rows = FALSE)
  } else {
    pheatmap::pheatmap(filled_data, angle_col = 45, border_color = "white",
                       color = hmColor, #breaks = myBreaks,
                       main = title)
  }
}


fill_symmetric_na <- function(mat) {
  is_na <- is.na(mat)
  mat[is_na] <- t(mat)[is_na]
  return(mat)
}

plot_type1_transformation <- function(fit_dat_sub, stats){
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
  return(filled_data)
}

plot_type2_transformation <- function(fit_dat, stats){
  fit_dat <- as.data.frame(fit_dat[,c("celltype1","celltype2","ROI",stats)])
  fit_dat[["pair"]] <- paste0(fit_dat[["celltype1"]],"-",fit_dat[["celltype2"]])
  fit_dat <- fit_dat[,c("pair","ROI",stats)]
  filled_data <- reshape(fit_dat, idvar = "pair", timevar = "ROI", 
                         direction = "wide") |>
    col2rownames("pair")
  colnames(filled_data) <- gsub(paste0("^",stats,"\\."),"ROI#",colnames(filled_data))
  return(filled_data)
}

