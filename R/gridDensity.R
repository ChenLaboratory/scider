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
#' @param isVisium Options of 'none','visium', and 'visiumHD'. If TRUE, converts
#' coordinates from pixel to um and fit the density grid to the same Visium spots 
#' arrangement. visium will use hexagonal while visiumHD will use rectangular grid.
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
                        id = if (isVisium!="none") NULL else "cell_type",
                        coi = NULL,
                        feature = NULL,
                        assay = "counts",
                        kernel = "gaussian",
                        bandwidth = NULL,
                        ngrid.x = NULL,
                        grid.length.x = NULL,
                        diggle = FALSE,
                        grid.type = c("hex", "square"),
                        isVisium = c("none","visium","visiumHD")
) {
  if(isTRUE(isVisium)) isVisium = "visium" # Backward compatibility
  else {isVisium <- match.arg(isVisium)}
  # Checks for Visium
  if (isVisium != "none") {
    if (is.null(spe$array_col) || 
        is.null(spe$array_row) ||
        is.null(spe$in_tissue)) {
      stop("Visium must have array_col, array_row, and in_tissue in colData")
    }
    if (isVisium == "visium") {
      if (missing(grid.type)) grid.type <- "hex"
      spe <- realignVisium(spe)
    } else { # visiumHD
      if (missing(grid.type)) grid.type <- "square"
      spe <- realignVisiumHD(spe)
      visium_bin_size <- .guessVisiumHDBin(spe)
      if (is.null(ngrid.x) && is.null(grid.length.x)) {
        grid.length.x <- 16 # visium_bin_size
      }
    }
  }
  
  grid.type <- match.arg(grid.type)
  
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
    w <- Matrix::t(SummarizedExperiment::assay(spe,assay)[feature,,drop=FALSE])
    weights <- cbind(weights,w)
  }
  # overall weight
  if (isVisium!="none") {
    weights <- cbind(weights,overall=spe$in_tissue)
  } else {
    weights <- cbind(weights,overall=rep.int(1,nrow(weights)))
  }
  clean_names <- paste("density",
                       janitor::make_clean_names(colnames(weights)),
                       sep="_")
  
  # Cells' coords
  spatialCoordsNames(spe)[1:2] <- c("x_centroid", "y_centroid")
  coord <- spatialCoords(spe)
  xlim <- range(coord[,"x_centroid"])
  ylim <- range(coord[,"y_centroid"])
  if (isVisium=="visium" && min(spe$array_col)%%2) {
    # Move xlim in case first columns is odd instead of even
    xlim[0] = xlim[0]-50
  } else if (isVisium=="visiumHD") {
    # Expand lims so that cells fall in the middle instead of the corner of bins
    xlim = xlim + visium_bin_size/2*c(-1,1)
    ylim = ylim + visium_bin_size/2*c(-1,1)
  }
  
  # Grid size. If both are provided, use grid.length.x
  if (!is.null(ngrid.x) && !is.null(grid.length.x)) {
    ngrid.x <- NULL
  }
  ngrid.x <- ngrid.x %||% (diff(xlim)/(grid.length.x %||% 100))
  
  # if (isVisium != "none") {
  # one_to_one <- FALSE
  # if (isVisium == "visium") {
  #   # n_col <- diff(range(spe$array_col))/`if`(isVisium=="visiumHD",1,2)
  #   n_col <- diff(range(spe$array_col))/2
  #   one_to_one <- isTRUE(all.equal(ngrid.x,n_col))
  #   if (!one_to_one) {
  #     #TODO: warning message not entirely accurate.
  #     message("For Visium, grid.length.x should be a divisible by 100 to exactly align each spot to a hexagon")
  #   }
  # }
  
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
  spe@metadata$grid_info$isVisium <- isVisium
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
  
  # if (filterToVisiumSpot && one_to_one) {
  #   hcellsInTissue <- hexDensity::xy2hcell(x=coord,
  #                                          xbins=spe@metadata$grid_info$xbins,
  #                                          xbnds=xlim,
  #                                          ybnds=ylim,
  #                                          shape=spe@metadata$grid_info$shape)
  #   hcellsInTissue <- sort(unique(hcellsInTissue)) #TODO: unique may not be needed
  #   spe@metadata$grid_density <- spe@metadata$grid_density[hcellsInTissue,]
  #   # sort cells to same order as gridpoints for subsetting
  #   spe=spe[,order(spe$array_row,spe$array_col)]
  #   spe@metadata$grid_info$gridLevelAnalysis <- TRUE
  # }
  
  return(spe)
}

