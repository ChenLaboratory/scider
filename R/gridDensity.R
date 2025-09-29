#' Perform kernel density estimation on SpatialExperiment for
#' cell types of interest
#'
#' @param spe A SpatialExperiment object.
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default. Set to NULL for overall density.
#' @param coi A character vector of cell types of interest (COIs).
#' Default to all cell types.
#' @param feature Feature(s) to calculate density with. Must be in rownames(spe).
#' @param assay Name of assay to use for finding feature(s).
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
                        id = if (isVisium) NULL else "cell_type",
                        coi = NULL,
                        feature = NULL,
                        assay = "counts",
                        kernel = "gaussian",
                        bandwidth = NULL,
                        ngrid.x = NULL,
                        grid.length.x = NULL,
                        diggle = FALSE,
                        grid.type = c("hex", "square"),
                        isVisium = FALSE,
                        filterToVisiumSpot = isVisium) {
  grid.type <- match.arg(grid.type)
  
  # Checks for Visium
  if (isVisium) {
    if (grid.type == "square") {
      grid.type <- "hex"
      message("Switching grid.type to hex for Visium")
    }
    if (is.null(spe$array_col) || 
        is.null(spe$array_row) ||
        is.null(spe$in_tissue)) {
      stop("Visium must have array_col, array_row, and in_tissue in colData")
    }
    
    # Scale the distances between points to 100 units and straighten the row/col
    spe <- realignVisium(spe)
  }
  
  weights <- matrix(ncol=0,nrow=ncol(spe))
  # id weight
  if (!is.null(id)) {
    if (!id %in% colnames(colData(spe))) { 
      message(paste(id, "is not a column of the colData. Skipping",id))
    }
    if (!all(coi %in% unique(spe@colData[[id]]))) {
      stop(paste(
        paste0(coi[!coi %in% unique(spe@colData[[id]])],collapse = ", "), 
        "not found in data!"))
    }
    
    if (is.numeric(spe[[id]])) {
      w <- matrix(spe[[id]],ncol=1,dimnames=list(NULL,id))
      w[is.na(w)] <- 0
    } else { # One-hot matrix
      f <- as.factor(spe[[id]])
      w <- matrix(0,nrow=ncol(spe),ncol=nlevels(f),dimnames=list(NULL,levels(f)))
      for (i in seq_along(f)) {
        w[i,f[i]] = 1
      }
      if (!is.null(coi)) w <- w[,coi,drop=FALSE]
    }
    weights <- cbind(weights,w)
  }
  # features weight
  f_not <- !(feature %in% rownames(spe))
  if (any(f_not)) {
    message(paste(paste0(feature[f_not],collapse = ", "),
                  "not found in rownames. Skipping them"))
    feature <- feature[!f_not]
  }
  if (!is.null(feature)) {
    w <- t(as.matrix(SummarizedExperiment::assay(spe,assay)[feature,,drop=FALSE]))
    # w <- t(as.matrix(spe@assays@data[[assay]][feature,,drop=FALSE]))
    weights <- cbind(weights,w)
  }
  # overall weight
  if (isVisium) {
    weights <- cbind(weights,overall=spe$in_tissue)
  } else {
    weights <- cbind(weights,overall=rep.int(1,nrow(weights)))
  }
  clean_names <- paste("density",
                       janitor::make_clean_names(colnames(weights)),
                       sep="_")
  
  # Cells' coords
  spatialCoordsNames(spe) <- c("x_centroid", "y_centroid")
  coord <- spatialCoords(spe)
  xlim <- range(coord[,"x_centroid"])
  ylim <- range(coord[,"y_centroid"])
  
  # Grid size. If both are provided, use grid.length.x
  if (!is.null(ngrid.x) && !is.null(grid.length.x)) {
    ngrid.x <- NULL
  }
  ngrid.x <- ngrid.x %||% (diff(xlim)/(grid.length.x %||% 100))
  
  
  if (isVisium) {
    one_to_one <- isTRUE(all.equal(ngrid.x,diff(range(spe$array_col))/2))
    if (!one_to_one) {
      message("For Visium, grid.length.x should be a divisible by 100 to exactly align each spot to a hexagon")
    }
  }
    
    
  
  # Calculate bandwidth
  if (is.null(bandwidth)) {
    if (!is.null(spe@metadata$grid_info$bandwidth)) {
      bandwidth <- spe@metadata$grid_info$bandwidth
      message("Reusing existing bandwidth for kernel smoothing!")
    } else {
      pts <- ppp(coord[, 1], coord[, 2], xlim, ylim)
      bandwidth <- bw.diggle(pts) * 4
    }
  }
  
  # Reset when the function is rerun again
  # spe@metadata$grid_density <- spe@metadata$grid_info <- NULL
  densFunc <- `if`(grid.type=="hex",computeDensityHex,computeDensity)
  # Set up info about the grid
  res <- densFunc(x = coord,
                  kernel = kernel,
                  bandwidth = bandwidth,
                  ngrid.x = ngrid.x,
                  xlim = xlim,
                  ylim = ylim,
                  gridInfo = TRUE)
  spe@metadata$grid_density <- res$grid_density
  spe@metadata$grid_info <- res$grid_info
  # Add densities
  for (ii in seq_len(ncol(weights))) {
    spe@metadata$grid_density <- cbind(
      spe@metadata$grid_density,
      densFunc(x = coord,
               kernel = kernel,
               bandwidth = bandwidth,
               weights = weights[,ii],
               ngrid.x = ngrid.x,
               xlim = xlim,
               ylim = ylim,
               diggle = diggle))
    colnames(spe@metadata$grid_density)[ncol(spe@metadata$grid_density)] = clean_names[[ii]]
  }
  
  if (grid.type=="hex") {
    # Filter grid_density to same as Visium spot.
    if (filterToVisiumSpot && isVisium && one_to_one) {
      hcellsInTissue <- hexDensity::xy2hcell(x=coord,
                                             xbins=spe@metadata$grid_info$xbins,
                                             xbnds=xlim,
                                             ybnds=ylim,
                                             shape=spe@metadata$grid_info$shape)
      hcellsInTissue <- sort(unique(hcellsInTissue))
      spe@metadata$grid_density <- spe@metadata$grid_density[hcellsInTissue,]
      spe@metadata$grid_info$gridLevelAnalysis <- TRUE
    }
    
    if (isVisium) spe@metadata$grid_info$isVisium <- TRUE
  }
  
  return(spe)
}

