#' Draw a contour region on some density level
#'
#' @param spe A SpatialExperiment object.
#' @param contour Name in metadata
#' @param coi A character vector of cell types of interest (COIs).
#' @param cutoff A numeric scalar specifying the density cutoff.
#'
#' @return An sf object of the contour region of the specified level. 
#' @export
#'
#' @examples
contour2sf <- function(spe, contour, coi, cutoff) {
  
  if(!requireNamespace("sf",quietly=TRUE)) stop("sf required but is not available")
  if(!requireNamespace("lwgeom",quietly=TRUE)) stop("lwgeom required but is not available")
  
  if(is.null(spe@metadata[[contour]])) stop("Have to calculate contours first!")
  
  clines <- as.data.frame(spe@metadata[[contour]])
  xlim <- c(min(clines$x), max(clines$x))
  ylim <- c(min(clines$y), max(clines$y))
  lims <- c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])
  canvas_sf <- sf::st_sf(sf::st_as_sfc(sf::st_bbox(lims)))
  levs <- sort(unique(clines$cutoff))
  lev_code <- findInterval(cutoff, levs, rightmost.closed = TRUE)
  clines_lev <- clines[clines$cutoff == cutoff, c("x", "y", "piece")]
  
  coi_clean <- janitor::make_clean_names(coi)
  dens_cols <- paste("density", coi_clean, sep="_")
  
  dens <- as.data.frame(spe@metadata$grid_density)
  dens$density_coi_average <- rowMeans(dens[, which(colnames(dens) %in% dens_cols), drop=FALSE])
  
  # grid centroid points
  dens$cutoff <- findInterval(dens$density_coi_average, levs, rightmost.closed = TRUE)
  grids_pts_sf <- dens[, c("node", "x_grid", "y_grid", "cutoff")]
  grids_pts_sf <- sf::st_as_sf(grids_pts_sf, coords = c("x_grid", "y_grid"))
  
  clines_lev_sf <- lapply(unique(clines_lev$piece), function(pp) {
    
    line_piece <- clines_lev[clines_lev$piece == pp, ]
    line_piece_sf <- sf::st_as_sf(line_piece, coords = c("x", "y")) |> sf::st_combine() |> 
      sf::st_multilinestring() |> sf::st_combine()
    
    # check if the region is at boundary
    bbox <- sf::st_bbox(line_piece_sf)
    bbox_sf <- sf::st_sf(sf::st_as_sfc(bbox))
    
    if (any(bbox == lims)) {
      regions <- lwgeom::st_split(bbox_sf, line_piece_sf) |> sf::st_collection_extract("POLYGON")
      inds <- sf::st_intersects(regions, grids_pts_sf)
      avglevel <- sapply(inds, function(ii) mean(grids_pts_sf$cutoff[ii]))
      area_up <- regions[avglevel >= lev_code, ]
      area_down <- regions[avglevel < lev_code, ]
      sf::st_geometry(area_up) <- "area"
      sf::st_geometry(area_down) <- "area"
      area_up <- sf::st_sf(area_up)
      area_down <- sf::st_sf(area_down)
      # bbox at boundaries
      bbox_boundary <- bbox_sf
    } else {
      area <- sf::st_cast(line_piece_sf, "POLYGON")
      area <- sf::st_sf(area)
      inds <- sf::st_intersects(area, grids_pts_sf)
      avglevel <- mean(grids_pts_sf$cutoff[unlist(inds)])
      if (avglevel >= lev_code) {
        area_up <- area; area_down <- NULL
      } else {
        area_up <- NULL; area_down <- area
      }
      bbox_boundary <- NULL
    }
    
    return(list(area_up = area_up, area_down = area_down, bbox_boundary = bbox_boundary))
    
  })
  
  areas_up <- do.call(rbind, lapply(clines_lev_sf, "[[", "area_up"))
  areas_down <- do.call(rbind, lapply(clines_lev_sf, "[[", "area_down"))
  bboxes_boundary <- do.call(rbind, lapply(clines_lev_sf, "[[", "bbox_boundary"))
  
  if (!is.null(areas_down)) {
    areas_up_union <- sf::st_union(areas_up)
    out1 <- sf::st_covered_by(areas_down, areas_up_union, sparse = FALSE)
    out2 <- sf::st_intersects(areas_down, areas_up_union, sparse = FALSE)
    out <- out1 | out2
    areas_down_out <- areas_down[c(out), ]
    ndown <- nrow(areas_down_out)
    areas <- areas_up_union
    for (dd in 1:ndown) {
      anyUp <- c(sf::st_intersects(areas_down_out[dd, ], sf::st_union(areas_up), sparse = FALSE) & 
                   !sf::st_touches(areas_down_out[dd, ], sf::st_union(areas_up), sparse = FALSE) &
                   !sf::st_covered_by(areas_down_out[dd, ], sf::st_union(areas_up), sparse = FALSE))
      # anyDownOverlap <- c(sf::st_intersects(areas_down_out[dd, ], areas_down_out, sparse = FALSE) & 
      #                       !sf::st_touches(areas_down_out[dd, ], areas_down_out, sparse = FALSE) & 
      #                       !sf::st_covered_by(areas_down_out[dd, ], areas_down_out, sparse = FALSE))
      # if (anyUp & sum(anyDownOverlap) > 0L) {
      if (anyUp) {
        areas_down_out_dd <- sf::st_difference(areas_down_out[dd, ], sf::st_union(areas_up))
      } else { areas_down_out_dd <- areas_down_out[dd, ] }
      areas <- sf::st_difference(areas, areas_down_out_dd)
    }
    # check if there is any missed area
    missed <- !sf::st_intersects(areas_up, areas, sparse = FALSE)
    if (any(missed)) {
      areas <- sf::st_union(sf::st_sf(areas), areas_up[missed, ])
    }
    # check if rest of the region still has up areas
    bboxes_boundary <- sf::st_union(bboxes_boundary)
    area_rest <- sf::st_difference(canvas_sf, bboxes_boundary)
    area_rest <- sf::st_cast(area_rest, "POLYGON")
    # boundary grids
    boundar_ind <- dens$x_grid == min(dens$x_grid) | dens$x_grid == max(dens$x_grid) | 
      dens$y_grid == min(dens$y_grid) | dens$y_grid == max(dens$y_grid)
    boundary_grids_pts_sf <- dens[boundar_ind, c("node", "x_grid", "y_grid", "cutoff")]
    boundary_grids_pts_sf <- sf::st_as_sf(boundary_grids_pts_sf, coords = c("x_grid", "y_grid"))
    inds <- sf::st_intersects(area_rest, boundary_grids_pts_sf, sparse = FALSE)
    inds <- sapply(1:nrow(inds), function (bb) {
      sum(boundary_grids_pts_sf$cutoff[inds[bb, ]] >= lev_code)
    })
    area_rest <- area_rest[inds > 0L, ]
    if (nrow(area_rest) > 0L) {
      area_rest <- sf::st_union(area_rest)
      anyIntersect_down <- st_intersects(areas_down_out, sparse = FALSE)
      anyIntersect_down <- anyIntersect_down[lower.tri(anyIntersect_down, diag = FALSE)]
      if (any(anyIntersect_down)) {
        ndown <- nrow(areas_down)
        for (dd in 1:ndown) {
          area_still_up <- sf::st_difference(area_rest, areas_down[dd, ])
        }
      } else {
        areas_down_combine <- sf::st_combine(areas_down)
        area_still_up <- sf::st_difference(area_rest, areas_down_combine)
      }
      areas <- sf::st_union(areas, area_still_up)
    }
    areas <- sf::st_union(areas)
  } else {
    areas <- sf::st_union(areas_up)
  }
  return(areas)
}


