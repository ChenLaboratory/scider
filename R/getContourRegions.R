#' Calculate areas between every two density levels
#'
#' @param spe A SpatialExperiment object.
#' @param contour_name Name of contour in spe@metadata
#'
#' @return A list of sf objects, each representing the region between
#' two contour density levels.
#'

getContourRegions <- function(spe, contour_name) {
  contour_data <- spe@metadata[[contour_name]]
  levs <- unique(contour_data$cutoff)
  nlevs <- length(levs)
  
  # calculate regions for each level
  if (is.null(spe@metadata[[contour_name]])) stop("Have to calculate
                                           contours first!")
  dens_cols <- S4Vectors::metadata(contour_data)$densities
  dens <- as.data.frame(spe@metadata$grid_density)
  dens$density_coi_average <- rowMeans(dens[, which(colnames(dens) %in%
                                                      dens_cols),
                                            drop = FALSE
  ])
  
  # Canvas.
  grid_info <- spe@metadata$grid_info
  if (grid_info$grid_type=="hex") {
    # Hexagonal canvas' sides are jagged.
    x_l <- rep_len(c(grid_info$xlim[1],grid_info$xlim[1]+grid_info$xstep/2),
                   grid_info$dims[2])
    x_r <- x_l + grid_info$xstep*(grid_info$dims[1]-1)
    y <- grid_info$ylim[1]+grid_info$ystep*(seq_len(grid_info$dims[2])-1)
    canvas_sf <- sf::st_polygon(list(matrix(c(
      x_l,rev(x_r),x_l[1],
      y,rev(y),y[1]
    ),ncol=2)))
  } else {# Square
    range_x <- spe@metadata$grid_info$xcol[c(1,length(spe@metadata$grid_info$xcol))]
    range_y <- spe@metadata$grid_info$yrow[c(1,length(spe@metadata$grid_info$yrow))]
    canvas_sf <- sf::st_polygon(list(matrix(c(
      range_x[c(1,1,2,2,1)],
      range_y[c(1,2,2,1,1)]),
      ncol=2)))
  }
  
  # Turn contour into sf. each multilinestring is a contour level
  pieces <- split.data.frame(as.matrix(contour_data[,c("x","y","piece")]),contour_data$level)
  pieces_sf <- sf::st_as_sfc(lapply(pieces, function(i){
    sf::st_multilinestring(split.data.frame(i[,c("x","y")],i[,"piece"]))
  }))

  # Split canvas into regions using all contours
  canvas_sf <- sf::st_snap(canvas_sf,pieces_sf,tolerance = spe@metadata$grid_info$xstep*0.01)
  regions <- lwgeom::st_split(canvas_sf,pieces_sf)|>
    sf::st_collection_extract("POLYGON")
  
  ## Assign regions to contour level.
  interval <- vector(length=length(pieces_sf))
  inds <- sf::st_intersects(regions, pieces_sf)
  n_touched <- lengths(inds)
  # Case 1: Regions intersect 2 contour levels
  interval[n_touched==2] <- sapply(inds[n_touched==2],min)
  # Case 2: Regions intersect only 1 contour levels
  grids_pts_sf <- sf::st_as_sf(dens, coords = c("x_grid","y_grid"))
  interval[n_touched==1] <- findInterval(
    sapply(sf::st_intersects(regions[n_touched==1], grids_pts_sf), 
           function(i) mean(grids_pts_sf$density_coi_average[i])
           ),
    levs)
  
  # Combine regions with the same contour level
  all_areas <- lapply(seq_len(nlevs), function(ii) {
    sf::st_as_sf(sf::st_combine(regions[interval==ii]))
  })
  names(all_areas) <- as.character(seq_len(nlevs))
  all_areas
}
