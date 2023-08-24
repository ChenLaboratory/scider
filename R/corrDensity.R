#' Model density correlation between two cell types
#'
#' @param spe A SpatialExperiment object.
#' @param celltype1 Cell type 1 to compare.
#' @param celltype2 Cell type 2 to compare.
#' @param cal_all Logical. If set to true, will model all cell types pairwise.
#' @param by_roi Logical. Plot facet by ROIs or not.
#' @param ngrid Integer. Threshold (minimum number of grids) used to filter small ROIs.
#' @param fit Character. Options are "spline" and "linear".
#' @param df Integer. Degrees of freedom of the spline fit. 
#' Default to 3 (i.e., a cubic spline fit).
#'
#' @return SpatialExperiment object.
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
corrDensity <- function(spe, celltype1 = NULL, celltype2 = NULL,
                        cal_all = FALSE, by_roi = TRUE, ngrid = 20, 
                        fit = c("spline","linear"), df = 3){
  
  if (!("grid_density" %in% names(spe@metadata))){
    stop("Please run gridDensity before using this function.")
  }
  
  if (!("roi" %in% names(spe@metadata))){
    stop("Please run findROI before using this function.")
  }
  
  if (length(fit) == 2){
    fit <- "spline"
  }
  
  dens_dat <- as.data.frame(spe@metadata$grid_density)
  rois <- as.data.frame(spe@metadata$roi)
  
  # filter rois
  kp <- which(table(rois$component) >= ngrid)
  rois_f <- rois[rois$component %in% kp, ]
  
  if (isTRUE(cal_all)){
    den_cols <- colnames(dens_dat)[grepl("density_", colnames(dens_dat))]
    
    model_results <- data.frame()
    for(i in seq(length(den_cols))){
      if (i == length(den_cols)){
        break
      }
      ct1 <- den_cols[i]
      for(j in seq(i+1, length(den_cols))){
        ct2 <- den_cols[j]
        model_result <- fit_model(dens_dat, rois_f, ct1, ct2, by_roi, fit)
        model_results <- rbind(model_results, model_result)
      }
    }
  } else {
    ct1 <- paste0("density_",janitor::make_clean_names(celltype1))
    ct2 <- paste0("density_",janitor::make_clean_names(celltype2))
    model_results <- fit_model(dens_dat, rois_f, ct1, ct2, by_roi, fit)
  }
  
  spe@metadata$model_result <- S4Vectors::DataFrame(model_results)
  return(spe)
}


fit_model <- function(dens_dat, rois_f, ct1, ct2, by_roi, fit){
  
  model_data_pre <- merge(rois_f, dens_dat, by.x = "members", 
                          by.y = "node", all.x = TRUE, sort = FALSE)
  model_data <- model_data_pre[,c("component", c(ct1, ct2))]
  
  colnames(model_data) <- c("component","x","y")
  
  cpnts <- unique(model_data$component)
  
  n1 <- janitor::make_clean_names(ct1, case = "sentence", 
                                  replace = c("density"= ""))
  n2 <- janitor::make_clean_names(ct2, case = "sentence", 
                                  replace = c("density"= ""))
  
  model_result <- data.frame("celltype1" = rep(n1, length(cpnts)),
                             "celltype2" = rep(n2, length(cpnts)),
                             "coefficient" = 0,
                             "correlation" = 0)
  if (isTRUE(by_roi)){
    if (fit == "spline"){
      for (i in seq(length(cpnts))){
        model <- lm(x ~ splines::ns(y, df = df), 
                    data = model_data[model_data$component == cpnts[i],])
        model_result <- get_model_result(model, model_result, i)
        model_result[i,"roi"] <- cpnts[i]
      }
    } else if (fit == "linear"){
      for (i in seq(length(cpnts))){
        model <- lm(y ~ x, 
                    data = model_data[model_data$component == cpnts[i],])
        model_result <- get_model_result(model, model_result, i, TRUE)
        model_result[i,"roi"] <- cpnts[i]
      }
    } 
  } else {
    model_result <- model_result[1,]
    if (fit == "spline"){
      model <- lm(x ~ splines::ns(y, df = df) + component, 
                  data = model_data)
      model_result <- get_model_result(model, model_result, 1)
    } else if (fit == "linear"){
      model <- lm(x ~ y + component,
                  data = model_data)
      model_result <- get_model_result(model, model_result, 1, TRUE)
    }
  }
  return(model_result)
}


get_model_result <- function(model, model_result, linenum, isLinear = FALSE){
  
  coeff <- as.numeric(coef(summary(model))[,"Estimate"][2])
  corr <- ifelse(coeff > 0, "positive", "negative")
  model_result[linenum,"coefficient"] <- coeff
  model_result[linenum,"correlation"] <- corr
  
  if (isTRUE(isLinear)){
    t_val <- as.numeric(coef(summary(model))[,"t value"][2])
    p_val <- as.numeric(coef(summary(model))[,"Pr(>|t|)"][2])
    model_result[linenum,"t"] <- t_val
    model_result[linenum,"p.val"] <- p_val
  }
  return(model_result)
}

