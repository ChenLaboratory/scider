#' Plot density correlation between two cell types
#'
#' @param spe A SpatialExperiment object.
#' @param celltype1 Cell type 1 to compare.
#' @param celltype2 Cell type 2 to compare.
#' @param by_roi Logical. Plot facet by ROIs or not.
#' @param threshold Integer. Threshold used to filter ROIs.
#' @param ... aesthetic mappings to pass to `ggplot2::aes()`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' 
#' data("xenium_bc_spe")
#' 
#' coi <- c("Breast cancer", "Fibroblasts")
#' 
#' spe <- gridDensity(spe, coi = coi)
#' 
#' spe <- findROI(spe, coi = coi, method = "walktrap")
#' 
#' plotDensCor(spe, celltype1 = "Breast cancer", celltype2 = "Fibroblasts")
#' 
plotDensCor <- function(spe, celltype1 = NULL, celltype2 = NULL,
                        by_roi = TRUE, threshold = 50, ...){
  
  rois <- metadata(spe)$components
  dens_dat <- metadata(spe)$grid_density
  
  # filter rois
  kp <- which(table(rois$component) > threshold)
  rois_f <- rois[rois$component %in% kp, ]
  
  # clean names
  ct1 <- paste0("density_",janitor::make_clean_names(celltype1))
  ct2 <- paste0("density_",janitor::make_clean_names(celltype2))
  
  plotdf <- left_join(rois_f, dens_dat, by = c("members"="node")) |>
    dplyr::select(c("component", all_of(c(ct1, ct2))))
  
  x <- rlang::sym(ct1)
  y <- rlang::sym(ct2)
  
  # extract aes
  aesmap <- rlang::enquos(...)
  # compute plot
  aesmap <- aesmap[!names(aesmap) %in% c("x", "y")] # remove x,y mappings if present
  
  # split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes <- vapply(aesmap, rlang::quo_is_symbolic, FUN.VALUE = logical(1))
    defaultmap <- lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap <- aesmap[is_aes]
  } else {
    defaultmap <- list()
  }
  
  # set some default plotting parameters
  if (is.null(defaultmap$shape)) {
    defaultmap$shape <- 21
  }
  
  if (is.null(defaultmap$alpha)) {
    defaultmap$alpha <- 0.8
  }
  
  if(isTRUE(by_roi)){
    
    if (is.null(defaultmap$fill)) {
      defaultmap$fill <- "royalblue"
    }
    
    p <- ggplot2::ggplot(plotdf, ggplot2::aes(!!x, !!y, !!!aesmap)) +
      do.call(ggplot2::geom_point, defaultmap) +
      geom_smooth(method='lm', formula = y ~ splines::ns(x, df = 3), 
                  color = "red", se = FALSE) +
      facet_wrap(~component, scales = "free", 
                 labeller = labeller(component = function(label) paste0("ROI #", label))) +
      theme_classic()
  } else {
    p <- ggplot2::ggplot(plotdf, ggplot2::aes(!!x, !!y, !!!aesmap)) +
      do.call(ggplot2::geom_point, defaultmap) +
      theme_classic()
  }
  
  return(p)
}




