#' Test for density correlation between two cell types.
#'
#' @param spe A SpatialExperiment object.
#' @param coi Character vector for cell types of interest for density 
#' correlation analysis. Default is NULL, which is to consider all cell types
#' previously calculated in the gridDensity() step. 
#' @param trace Logical. If TRUE (default), print process pf testing. 
#'
#' @return A DataFrame containing the testing results.
#' @export
#' 
#' @importFrom SpatialPack modified.ttest
#' @importFrom stats pt
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
#' result <- corDensity(spe)
#'
corDensity <- function(spe, coi = NULL, trace = TRUE) {
  if (!("grid_density" %in% names(spe@metadata))) {
    stop("Please run gridDensity before using this function.")
  }
  
  if (!("roi" %in% names(spe@metadata))) {
    stop("Please run findROI before using this function.")
  }
  
  dens_dat <- as.data.frame(spe@metadata$grid_density)
  rois <- as.data.frame(spe@metadata$roi)
  
  # get cell type info
  den_cols <- colnames(dens_dat)[grepl("density_", colnames(dens_dat))]
  
  if(any(!is.null(coi))){
    coi_clean <- paste0("density_", janitor::make_clean_names(coi))
    coi.exist <- coi_clean %in% den_cols
    if(any(!coi.exist)) 
      stop(paste0(paste(coi[!coi.exist], collapse=", "), " not found in the data."))
    den_cols <- den_cols[den_cols %in% coi_clean]
  }
  
  nCT <- length(den_cols)
  if (nCT == 1) stop("Only one cell type detected.")
  
  # construct data table
  model_data <- merge(rois, dens_dat,
                      by.x = "members",
                      by.y = "node", all.x = TRUE, sort = FALSE
  )
  component <- factor(model_data$component)
  nGrids <- table(component)
  cpnts <- names(nGrids)
  nCpnts <- length(cpnts)
  
  nRows <- choose(nCT, 2) * nCpnts
  result.ROI <- data.frame(
    "celltype1" = rep("", nRows),
    "celltype2" = rep("", nRows),
    "ROI" = 0,
    "ngrid" = 0,
    "cor.coef" = 0,
    "t" = 0,
    "df" = 0,
    "p.Pos" = 0,
    "p.Neg" = 0
  )
  
  result.overall <- result.ROI[seq_len(choose(nCT, 2)), c(1, 2, 5, 8, 9)]
  
  for (i in seq_len((nCT - 1))) {
    for (j in (i + 1):nCT) {
      ct1 <- den_cols[i]
      ct2 <- den_cols[j]
      n1 <- janitor::make_clean_names(ct1,
                                      case = "sentence",
                                      replace = c("density" = "")
      )
      n2 <- janitor::make_clean_names(ct2,
                                      case = "sentence",
                                      replace = c("density" = "")
      )
      m <- choose(nCT, 2) - (choose(nCT - j, 1) + choose(nCT - i, 2))
      
      for (k in seq_len(nCpnts)) {
        if (trace) cat(paste("i =", i, ", j =", j, ", ROI", k, "\n"))
        data <- model_data[model_data$component == cpnts[k], c(ct1, ct2, "x", "y")]
        
        res <- modified.ttest(x=data[,1], y=data[,2], 
                              coords=data[,c("x","y")], nclass=7)
        tstat <- sqrt( res$Fstat * res$dof ) * sign(res$corr)
        
        ind <- (m - 1) * nCpnts + k
        result.ROI[ind, "celltype1"] <- n1
        result.ROI[ind, "celltype2"] <- n2
        result.ROI[ind, "ROI"] <- cpnts[k]
        result.ROI[ind, "ngrid"] <- nGrids[k]
        result.ROI[ind, "cor.coef"] <- res$corr
        result.ROI[ind, "t"] <- tstat
        result.ROI[ind, "df"] <- res$dof
        result.ROI[ind, "p.Pos"] <- pt(tstat, df=res$dof, lower.tail = FALSE)
        result.ROI[ind, "p.Neg"] <- pt(tstat, df=res$dof, lower.tail = TRUE)
      }
      
      # summarize test results across ROIs
      sel <- (m - 1) * nCpnts + seq_len(nCpnts)
      out.cor.coef <- result.ROI$cor.coef[sel]
      out.p.Pos <- result.ROI$p.Pos[sel]
      out.p.Neg <- result.ROI$p.Neg[sel]
      
      result.overall[m, "celltype1"] <- n1
      result.overall[m, "celltype2"] <- n2
      result.overall[m, "cor.coef"] <- tanh(mean(atanh(out.cor.coef)))
      result.overall[m, "p.Pos"] <- pchisq(-2 * sum(log(out.p.Pos)),
                                           df = 2 * nCpnts, lower.tail = FALSE)
      result.overall[m, "p.Neg"] <- pchisq(-2 * sum(log(out.p.Neg)),
                                           df = 2 * nCpnts, lower.tail = FALSE)
    }
  }
  
  result.ROI <- S4Vectors::DataFrame(result.ROI)
  result.overall <- S4Vectors::DataFrame(result.overall)
  
  return(list(ROI = result.ROI, overall = result.overall))
}
