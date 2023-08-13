#' Draw a contour region on some density level
#'
#' @param spe A SpatialExperiment object.
#' @param contour Name in metadata
#' @param coi A character vector of cell types of interest (COIs).
#' @param level 
#'
#' @return An sf object of the contour region of the specified level. 
#' @export
#'
#' @examples

contour2sf <- function(spe, contour, coi, level) {
  
  if (is.null(spe@metadata[[contour]]))
    stop("Have to calculate contours first!")
  
  clines <- spe@metadata[[contour]]
  xlim <- c(min(clines$x), max(clines$x))
  ylim <- c(min(clines$y), max(clines$y))
  lims <- c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])
  canvas_sf <- st_sf(st_as_sfc(st_bbox(lims)))
  levs <- sort(unique(clines$level))
  lev_code <- findInterval(level, levs, rightmost.closed = TRUE)
  clines_lev <- clines[clines$level == level, c("x", "y", "piece")]
  
  coi_clean <- janitor::make_clean_names(coi)
  dens_cols <- paste("density", coi_clean, sep="_")
  
  dens <- spe@metadata$grid_density
  dens$density_coi_average <- rowMeans(dens[, which(colnames(dens) %in% dens_cols), drop=FALSE])
  
  # grid centroid points
  dens$level <- findInterval(dens$density_coi_average, levs, rightmost.closed = TRUE)
  grids_pts_sf <- dens[, c("node", "x_grid", "y_grid", "level")]
  grids_pts_sf <- st_as_sf(grids_pts_sf, coords = c("x_grid", "y_grid"))
  
  clines_lev_sf <- lapply(unique(clines_lev$piece), function(pp) {
    
    line_piece <- clines_lev[clines_lev$piece == pp, ]
    line_piece_sf <- st_as_sf(line_piece, coords = c("x", "y")) |> st_combine() |> 
      st_multilinestring() |> st_combine()
    
    # check if the region is at boundary
    bbox <- st_bbox(line_piece_sf)
    bbox_sf <- st_sf(st_as_sfc(bbox))
    if (any(bbox == lims)) {
      regions <- st_split(bbox_sf, line_piece_sf) |> st_collection_extract("POLYGON")
      inds <- st_intersects(regions, grids_pts_sf)
      avglevel <- sapply(inds, function(ii) mean(grids_pts_sf$level[ii]))
      area_up <- regions[avglevel >= lev_code, ]
      area_down <- regions[avglevel < lev_code, ]
      st_geometry(area_up) <- "area"
      st_geometry(area_down) <- "area"
      area_up <- st_sf(area_up)
      area_down <- st_sf(area_down)
    } else {
      area <- st_cast(line_piece_sf, "POLYGON")
      area <- st_sf(area)
      inds <- st_intersects(area, grids_pts_sf)
      avglevel <- mean(grids_pts_sf$level[unlist(inds)])
      if (avglevel >= lev_code) {
        area_up <- area; area_down <- NULL
      } else {
        area_up <- NULL; area_down <- area
      }
    }
    
    return(list(area_up = area_up, area_down = area_down))
    
  })
  
  areas_up <- do.call(rbind, lapply(clines_lev_sf, "[[", "area_up"))
  areas_down <- do.call(rbind, lapply(clines_lev_sf, "[[", "area_down"))
  
  if (!is.null(areas_down)) {
    areas_up_union <- st_union(areas_up)
    out <- st_covered_by(areas_down, areas_up_union, sparse = FALSE)
    areas_down <- areas_down[c(out), ]
    areas_down <- st_combine(areas_down)
    areas <- st_difference(areas_up_union, areas_down)
    # check if there is any missed area
    missed <- !st_intersects(areas_up, areas, sparse = FALSE)
    if (any(missed)) {
      areas <- st_union(st_sf(areas), areas_up[missed, ])
    }
  } else {
    areas <- st_union(areas_up)
  }
  
  return(areas)
}


