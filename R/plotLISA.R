#' Plotting LISA (e.g. moran)
#' @param spe A SpatialExperiment object.
#' @param lisa Output from localMoran
#' @param overlay Option of grid or point. Depend on whether 
#' \link[scider]{localMoran} is calculated at grid or point.
#' @param type Option of cluster or logpvalue for plotting lisa's cluster or 
#' p-value, respectively.
#' @param ... Parameters pass to \link[scider]{plotGrid} or 
#' \link[scider]{plotSpatial}, depending on overlay.
#' 
#' @return a ggplot object
#' @examples 
#' data("xenium_bc_spe")
#' spe <- gridDensity(spe, coi = "Breast cancer")
#' dat <- spe@metadata$grid_density$density_breast_cancer
#' spe <- findNbrsGrid(spe)
#' res <- localMoran(spe,data1=dat, at = "grid")
#' plotLISA(spe,lisa = res)
#' 
#' @export
plotLISA <- function(spe, lisa,
                     overlay=c("grid","point"),
                     type=c("cluster","logpvalue"),
                     ...) {
  # Checking if overlay is valid
  n_point <- nrow(spe@colData)
  n_grid <- nrow(spe@metadata$grid_density) %||% -1
  if (missing(overlay)) {
    if (length(lisa$lisa)==n_point) overlay<-"point"
    else if (length(lisa$lisa)==n_grid) overlay<-"grid"
    else stop("Length of lisa must be equal to either number of cells or grid points")
  } else {
    overlay <- match.arg(overlay)
    if (overlay=="grid" && length(lisa$lisa)!=n_grid && 
        length(lisa$lisa)==n_point) {
      overlay <- "point"
      message("Length of lisa not equal to number of grid points. Switching
              overlay to 'point'")
    } else if (overlay=="point" && length(lisa$lisa)!=n_point && 
               length(lisa$lisa)==n_grid) {
      overlay <- "grid"
      message("Length of lisa not equal to number of cells Switching
              overlay to 'grid'")
    } else {
      stop("Length of lisa must be equal to either number of cells or grid points")
    }
  }

  keep <- lisa$cluster!="Undefined"
  
  # Getting relevant data from lisa 
  type <- match.arg(type)
  if (type=="cluster") {
    dat <- lisa$cluster[keep]
    cols <- col.lisa[sort(unique(dat))]
    cols.scale <- NULL
  } else {
    dat <- log(lisa$p.value[keep])
    cols <- col.pval
    cols.scale <- c(0,1-log(0.05)/min(dat),1)
  }
  
  # Plot
  if (overlay=="grid") {
    spe@metadata$grid_density <- spe@metadata$grid_density[keep,]
    spe@metadata$grid_density$dat <- dat
    p <- plotGrid(spe,group.by="dat",cols=cols,label=type,
                  cols.scale = cols.scale, ...)
  } else {
    spe <- spe[,keep]
    spe$dat <- dat
    p <- plotSpatial(spe,group.by="dat",cols=cols,label=type,
                     cols.scale = cols.scale, ...)
    p <- update_bound(p,
                      x=SpatialExperiment::spatialCoords(spe)[,1],
                      y=SpatialExperiment::spatialCoords(spe)[,2])
  }
  return(p)
}
