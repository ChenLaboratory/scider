#' Perform kernel density estimation on SpatialExperiment for
#' cell types of interest
#'
#' @param spe A SpatialExperiment object.
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default.
#' @param coi A character vector of cell types of interest (COIs).
#' Default to all cell types.
#' @param kernel The smoothing kernel. Options are "gaussian",
#' "epanechnikov", "quartic" or "disc". For hexagonal grid, only Gaussian is implemented
#' @param bandwidth The smoothing bandwidth. By default performing
#' automatic bandwidth selection using cross-validation using
#' function spatstat.explore::bw.diggle.
#' @param ngrid.x Number of grids in the x-direction. Ignored when
#' 'grid.length.x' is specified. Default to NULL.
#' @param grid.length.x Grid length in the x-direction. If both 
#' 'ngrid.x' and 'grid.length.x' are NULL, then 'grid.length.x'
#' is set to 100 (micron) by default.
#' @param diggle Logical. If TRUE, use the Jones-Diggle improved edge
#' correction. See spatstat.explore::density.ppp() for details.
#' @param grid.type Type of grid can be either hexagon or square.
#' @param isVisium Logical. If TRUE, fit hexagonal grids to Visium spots by 
#' replacing spatial coords with array rows & array cols. 
#' @param filterToVisiumSpot Logical. If TRUE, filter grid polygons to only 
#' those with a Visium spot underneath.
#'
#' @return A SpatialExperiment object. Grid density estimates for
#' all cell type of interest are stored in spe@metadata$grid_density.
#' Grid information is stored in spe@metadata$grid_info
#'
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe <- gridDensity(spe)
#'
gridDensity <- function(spe,
                        id = if (isVisium) "in_tissue" else "cell_type",
                        coi = NULL,
                        kernel = "gaussian",
                        bandwidth = NULL,
                        ngrid.x = NULL,
                        grid.length.x = NULL,
                        diggle = FALSE,
                        grid.type = c("hex", "square"),
                        isVisium = FALSE,
                        filterToVisiumSpot = isVisium) {
  grid.type <- match.arg(grid.type)
  if (isVisium && grid.type == "square") {
    grid.type <- "hex"
    message("Switching grid.type to hex for Visium")
  }
  
  if (!id %in% colnames(colData(spe))) {
    stop(paste(id, "is not a column of the colData."))
  }
  
  # (Default COI for non-visium) OR (visium without default id)
  if ((is.null(coi) && !isVisium) || (isVisium && id != "in_tissue")) {
    coi <- names(table(colData(spe)[[id]]))
  }
  
  if (length(which(!coi %in% names(table(colData(spe)[[id]])))) > 0L) {
    stop(paste(paste0(
      coi[which(!coi %in%
                  names(table(colData(spe)[[id]])))],
      collapse = ", "
    ), "not found in data!", sep = " "))
  }

  coi <- c(coi, "overall")
  coi_clean <- janitor::make_clean_names(coi)
  
  # define canvas
  if (isVisium) {

    if (is.null(colData(spe)$array_col) || 
        is.null(colData(spe)$array_row) ||
        is.null(colData(spe)$in_tissue)) {
      stop("Visium must have array_col, array_row, and in_tissue in colData")
    }
    spatialCoords(spe) <- cbind((colData(spe)$array_col)*50,
                                (colData(spe)$array_row)*50*sqrt(3))
  }
  spatialCoordsNames(spe) <- c("x_centroid", "y_centroid")
  coord <- spatialCoords(spe)
  xlim <- range(coord[,"x_centroid"])
  ylim <- range(coord[,"y_centroid"])
  
  # Calculate bandwidth
  pts <- ppp(coord[, 1], coord[, 2], xlim, ylim)
  if (is.null(bandwidth) & !is.null(spe@metadata$grid_info$bandwidth)) {
    bandwidth <- spe@metadata$grid_info$bandwidth
    message("Reusing existing bandwidth for kernel smoothing!")
  }
  if (is.null(bandwidth)) {
    bandwidth <- bw.diggle(pts) * 4
  }
  
  if (is.null(spe@metadata)) spe@metadata <- list()

  if(is.null(ngrid.x) && is.null(grid.length.x)) 
    grid.length.x <- 100

  # Reset when the function is rerun again
  spe@metadata$grid_density <- spe@metadata$grid_info <- NULL
  
  # compute density for each cell type and then, filter
  if (grid.type=="hex") {
    for (ii in seq_len(length(coi))) {
        if(coi[ii] != "overall"){
            # subset data to this COI
            sub <- which(colData(spe)[[id]] == coi[ii])
            obj <- spe[, sub]
        } else if (isVisium) {
          obj <- spe[, which(colData(spe)$in_tissue == 1)]
        } else 
            obj <- spe

      # compute density
      out <- computeDensityHex(obj,
                            kernel = kernel,
                            bandwidth = bandwidth,
                            ngrid.x = ngrid.x,
                            grid.length.x = grid.length.x,
                            xlim = xlim, ylim = ylim, diggle = diggle,
                            isVisium = isVisium
      )
      RES <- out$grid_density
      
      if (is.null(spe@metadata$grid_density)) {
        spe@metadata <- list("grid_density" = RES[, seq_len(4)])
        spe@metadata$grid_density$node <- paste(
          spe@metadata$grid_density$node_x,
          spe@metadata$grid_density$node_y,
          sep = "-"
        )
      }
  
      spe@metadata$grid_density <- cbind(spe@metadata$grid_density,RES$density)
      colnames(spe@metadata$grid_density)[5 + ii] <- paste("density",
                                                           coi_clean[ii],
                                                           sep = "_"
      )
      
      # grid info
      if (is.null(spe@metadata$grid_info)) {
        spe@metadata$grid_info <- list(
          dims = out$density_est@dimen[2:1],
          xlim=xlim,
          ylim=ylim,
          xstep=diff(xlim)/out$density_est@xbins,
          ystep=(diff(ylim)*sqrt(3))/(2*out$density_est@shape*out$density_est@xbins),
          xbins=out$density_est@xbins,
          shape=out$density_est@shape,
          bandwidth=bandwidth,
          grid_type = "hex"
        )
      }
    }
    #Filter grid_density to same as Visium spot.
    if (filterToVisiumSpot==TRUE && isVisium == TRUE && grid.length.x==100) {
      hcellsInTissue <- hexDensity::xy2hcell(x=spatialCoords(spe),
                                             xbins=out$density_est@xbins,
                                             xbnds=xlim,
                                             ybnds=ylim,
                                             shape=out$density_est@shape)
      spe@metadata$grid_density <- spe@metadata$grid_density[hcellsInTissue,]
      spe@metadata$grid_info$gridLevelAnalysis <- TRUE
    }
    if (isVisium==TRUE) spe@metadata$grid_info$isVisium <- TRUE
  } else {
    for (ii in seq_len(length(coi))) {

        if(coi[ii] != "overall"){
            # subset data to this COI
            sub <- which(colData(spe)[[id]] == coi[ii])
            obj <- spe[, sub]
        } else 
            obj <- spe

      # compute density
      out <- computeDensity(obj,
                               mode = "pixels", kernel = kernel,
                               bandwidth = bandwidth,
                               ngrid.x = ngrid.x,
                               grid.length.x = grid.length.x,
                               xlim = xlim, ylim = ylim, diggle = diggle
      )
      RES <- out$grid_density
      
      ngrid.x <- out$density_est$dim[2]
      ngrid.y <- out$density_est$dim[1]
      
      if (is.null(spe@metadata$grid_density)) {
        spe@metadata <- list("grid_density" = RES[, seq_len(2)])
        # horizontal ind
        spe@metadata$grid_density$node_x <- rep(seq_len(ngrid.x),
                                                each = ngrid.y
        )
        # vertical ind
        spe@metadata$grid_density$node_y <- rep(
          seq_len(ngrid.y),
          ngrid.x
        )
        spe@metadata$grid_density$node <- paste(
          spe@metadata$grid_density$node_x,
          spe@metadata$grid_density$node_y,
          sep = "-"
        )
      }
      
      spe@metadata$grid_density <- cbind(spe@metadata$grid_density,RES$density)
      colnames(spe@metadata$grid_density)[5 + ii] <- paste("density",
                                                           coi_clean[ii],
                                                           sep = "_"
      )
      
      # grid info
      if (is.null(spe@metadata$grid_info)) {
        spe@metadata$grid_info <- list(
          dims = c(ngrid.x, ngrid.y),
          xlim = xlim,
          ylim = ylim,
          xcol = out$density_est$xcol,
          yrow = out$density_est$yrow,
          xstep = out$density_est$xstep,
          ystep = out$density_est$ystep,
          bandwidth = bandwidth,
          grid_type = "square"
        )
      }
    }
  }

  return(spe)
}
