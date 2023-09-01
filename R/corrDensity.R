#' Test for density correlation between two cell types.
#'
#' @param spe A SpatialExperiment object.
#' @param ngrid Integer. Threshold (minimum number of grids) used to filter small ROIs.
#' @param by.roi Logical. If TRUE (default), then return the testing results at ROI level. 
#' If FALSE, then combine the testing results across all ROIs.
#'
#' @return A DataFrame containing the testing results.
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
#' result <- corrDensity(spe, ngrid = 30)
#' 
corrDensity <- function(spe, ngrid = 20, by.roi = TRUE){

  if (!("grid_density" %in% names(spe@metadata))){
    stop("Please run gridDensity before using this function.")
  }
  
  if (!("roi" %in% names(spe@metadata))){
    stop("Please run findROI before using this function.")
  }

  dens_dat <- as.data.frame(spe@metadata$grid_density)
  rois <- as.data.frame(spe@metadata$roi)

  # get cell type info
  den_cols <- colnames(dens_dat)[grepl("density_", colnames(dens_dat))]
  nCT <- length(den_cols)
  if (nCT == 1) stop("Only one cell type detected.")
  
  # filter rois
  kp <- which(table(rois$component) >= ngrid)
  rois_f <- rois[rois$component %in% kp, ]

  # construct data table
  model_data <- merge(rois_f, dens_dat, by.x = "members", 
                          by.y = "node", all.x = TRUE, sort = FALSE)
  component <- factor(model_data$component)
  nGrids <- table(component)
  cpnts <- names(nGrids)
  nCpnts <- length(cpnts)

  nRows <- choose(nCT,2)*nCpnts
  result.ROI <- data.frame("celltype1" = rep("", nRows),
                       "celltype2" = rep("", nRows),
                       "ROI" = 0,
                       "ngrid" = 0,
                       "cor.coef" = 0,
                       "t" = 0, 
                       "df" = 0,
                       "p.Pos" = 0,
                       "p.Neg" = 0)

  result.overall <- result.ROI[1:choose(nCT,2), c(1,2,5,8,9)]

  for(i in 1:(nCT-1)) for(j in (i+1):nCT){
    ct1 <- den_cols[i]
    ct2 <- den_cols[j]
    n1 <- janitor::make_clean_names(ct1, case = "sentence", replace = c("density"= ""))
    n2 <- janitor::make_clean_names(ct2, case = "sentence", replace = c("density"= ""))
    m <- choose(nCT, 2) - (choose(nCT-j, 1) + choose(nCT-i, 2))

    for(k in 1:nCpnts){
      res.pos <- cor.test(~ get(ct1) + get(ct2), data = model_data, 
        alternative = "greater", subset = model_data$component == cpnts[k])
      res.neg <- cor.test(~ get(ct1) + get(ct2), data = model_data, 
        alternative = "less", subset = model_data$component == cpnts[k])
      ind <- (m-1)*nCpnts + k

      result.ROI[ind, "celltype1"] <- n1
      result.ROI[ind, "celltype2"] <- n2
      result.ROI[ind, "ROI"] <- cpnts[k]
      result.ROI[ind, "ngrid"] <- nGrids[k]
      result.ROI[ind, "cor.coef"] <- res.pos$estimate
      result.ROI[ind, "t"] <- res.pos$statistic
      result.ROI[ind, "df"] <- res.pos$parameter
      result.ROI[ind, "p.Pos"] <- res.pos$p.value
      result.ROI[ind, "p.Neg"] <- res.neg$p.value
    }

    # summarize test results across ROIs    
    sel <- (m-1)*nCpnts + 1:nCpnts
    out.cor.coef <- result.ROI$cor.coef[sel] 
    out.p.Pos <- result.ROI$p.Pos[sel] 
    out.p.Neg <- result.ROI$p.Neg[sel] 

    result.overall[m, "celltype1"] <- n1
    result.overall[m, "celltype2"] <- n2
    result.overall[m, "cor.coef"] <- tanh(mean(atanh(out.cor.coef)))
    result.overall[m, "p.Pos"] <- pchisq(-2*sum(log(out.p.Pos)), df=2*nCpnts, lower.tail = FALSE)
    result.overall[m, "p.Neg"] <- pchisq(-2*sum(log(out.p.Neg)), df=2*nCpnts, lower.tail = FALSE)
  }

  result.ROI <- S4Vectors::DataFrame(result.ROI)
  result.overall <- S4Vectors::DataFrame(result.overall)

  if(by.roi)
    return(result.ROI)
  else
    return(result.overall)
}
