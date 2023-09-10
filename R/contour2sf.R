#' Draw a contour region on some density level
#'
#' @param spe A SpatialExperiment object.
#' @param contour Name in metadata.
#' @param coi A character vector of cell types of interest (COIs).
#' @param cutoff A numeric scalar specifying the density cutoff.
#'
#' @return An sf object of the contour region of the specified level. 
#' 
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
  grids_pts_sf <- dens[, c("node", "x_grid", "y_grid", "density_coi_average")]
  grids_pts_sf <- sf::st_as_sf(grids_pts_sf, coords = c("x_grid", "y_grid"))
  
  # all pieces
  all_clines_sf <- lapply(unique(clines_lev$piece), function(pp) {
    line_piece <- clines_lev[clines_lev$piece == pp, ]
    sf::st_as_sf(line_piece, coords = c("x", "y")) |> sf::st_combine() |> 
      sf::st_multilinestring() |> sf::st_combine()
  })
  names(all_clines_sf) <- unique(clines_lev$piece)
  
  clines_lev_sf <- lapply(names(all_clines_sf), function(pp) {
    
    line_piece_sf <- all_clines_sf[[pp]]
    # check if the region is at boundary
    bbox <- sf::st_bbox(line_piece_sf)
    bbox_sf <- sf::st_sf(sf::st_as_sfc(bbox))
    
    if (any(bbox == lims)) {
      # if the line intersects with parallel boundaries
      touch_which_boundary <- substr(names(bbox)[bbox == lims], start = 1, stop = 1)
      if (length(touch_which_boundary) >= 2L & any(table(touch_which_boundary) == 2L)) {
        bbox_parabound <- bbox_sf
      } else {
        bbox_parabound <- NULL
      }
      regions <- lwgeom::st_split(bbox_sf, line_piece_sf) |> sf::st_collection_extract("POLYGON")
      other_clines_ind <- which(names(all_clines_sf) != as.character(pp))
      if (length(other_clines_ind) > 0L) {
        other_clines_sf <- do.call(rbind, lapply(all_clines_sf[other_clines_ind], sf::st_sf))
        whether_cross <- c(sf::st_crosses(other_clines_sf, bbox_sf, sparse = FALSE))
        if (any(whether_cross)) {
          sub_region_ind <- sf::st_intersects(regions, other_clines_sf[whether_cross, ], 
                                              sparse = FALSE)
          sub_region_ind <- which(rowSums(sub_region_ind) > 0L)
          subregions <- lapply(sub_region_ind, function(subr) {
            lwgeom::st_split(regions[sub_region_ind, ], 
                             sf::st_combine(other_clines_sf[whether_cross, ])) |>
              sf::st_collection_extract("POLYGON")
          })
          regions <- regions[-sub_region_ind, ]
          regions <- rbind(do.call(rbind, subregions), regions)
        }
      } 
      inds <- sf::st_intersects(regions, grids_pts_sf)
      avglevel <- sapply(inds, function(ii) mean(grids_pts_sf$density_coi_average[ii]))
      avglevel <- findInterval(avglevel, levs, rightmost.closed = TRUE)
      area_up <- regions[avglevel >= lev_code, ]
      area_down <- regions[avglevel < lev_code, ]
      sf::st_geometry(area_up) <- "area"
      sf::st_geometry(area_down) <- "area"
      area_up <- sf::st_sf(area_up)
      area_down <- sf::st_sf(area_down)
      cline_touches_bound <- TRUE
    } else {
      area <- sf::st_cast(line_piece_sf, "POLYGON")
      area <- sf::st_sf(area)
      inds <- sf::st_intersects(area, grids_pts_sf)
      avglevel <- mean(grids_pts_sf$density_coi_average[unlist(inds)])
      avglevel <- findInterval(avglevel, levs, rightmost.closed = TRUE)
      if (avglevel >= lev_code) {
        area_up <- area; area_down <- NULL; bbox_parabound <- NULL
      } else {
        area_up <- NULL; area_down <- area; bbox_parabound <- NULL
      }
      cline_touches_bound <- FALSE
    }
    
    if (!is.null(area_up)) {
      area_up <- area_up[!is.na(sf::st_dimension(area_up)), ]
      if (nrow(area_up) == 0L) area_up <- NULL
    }
    if (!is.null(area_down)) {
      area_down <- area_down[!is.na(st_dimension(area_down)), ]
      if (nrow(area_down) == 0L) area_down <- NULL
    }
    
    return(list(area_up = area_up, 
                area_down = area_down, 
                cline_touches_bound = cline_touches_bound, 
                bbox_parabound = bbox_parabound))
    
  })
  
  areas_up <- do.call(rbind, lapply(clines_lev_sf, "[[", "area_up"))
  areas_down <- do.call(rbind, lapply(clines_lev_sf, "[[", "area_down"))
  cline_touches_bound <- unlist(lapply(clines_lev_sf, "[[", "cline_touches_bound"))
  bbox_parabound <- do.call(rbind, lapply(clines_lev_sf, "[[", "bbox_parabound"))
  
  if (!is.null(areas_down)) {
    areas_up_union <- sf::st_union(areas_up)
    out <- sf::st_covered_by(areas_down, areas_up_union, sparse = FALSE)
    areas_down_out <- areas_down[c(out), ]
    if (nrow(areas_down_out) > 0L) {
      # flatten out overlaps
      downHasOverlap <- sf::st_intersection(areas_down_out)
      areas_down_out <- sf::st_sf(sf::st_geometry(downHasOverlap))
      # remove MULTILINESTRING
      areas_down_out <- areas_down_out[sf::st_geometry_type(areas_down_out) %in% c("POLYGON",
                                                                                   "MULTIPOLYGON"), ]
      # check if downs are really down
      areas_down_out_code <- sapply(1:nrow(areas_down_out), function(xx) {
        xx <- areas_down_out[xx, ]
        inds <- sf::st_intersects(xx, grids_pts_sf)
        avglevel <- sapply(inds, function(ii) mean(grids_pts_sf$density_coi_average[ii]))
        avglevel <- findInterval(avglevel, levs, rightmost.closed = TRUE)
        return(avglevel)
      })
      areas_down_out <- areas_down_out[areas_down_out_code < lev_code, ]
      areas_down_out <- sf::st_combine(areas_down_out)
      areas <- sf::st_difference(areas_up_union, areas_down_out)
      # check if there is any missed area
      missed_up <- !sf::st_intersects(areas_up, areas, sparse = FALSE)
      if (any(missed_up)) {
        areas <- sf::st_union(sf::st_sf(areas), areas_up[missed_up, ])
      }
    } else {
      areas <- sf::st_union(areas_up)
    }
  } else {
    areas <- sf::st_union(areas_up)
  }
  
  # if there are clines cross parallel boundaries of the canvas
  if (!is.null(bbox_parabound)) {
    # get regions outside the bbox(es)
    bbox_parabound <- sf::st_union(bbox_parabound)
    stripes <- sf::st_difference(canvas_sf, bbox_parabound)
    if (nrow(stripes) > 0L) {
      stripes <- sf::st_cast(stripes, to = "POLYGON")
      # for each polygon, remove the existing up areas, and see if the rest of the region still has up regions
      stripes_up <- lapply(1:nrow(stripes), function(ss) {
        this_stripe <- stripes[ss, ]
        this_stripe_minus_up <- sf::st_difference(this_stripe, areas)
        this_stripe_minus_up <- sf::st_union(this_stripe_minus_up)
        this_stripe_minus_up <- sf::st_cast(this_stripe_minus_up, to = "POLYGON")
        this_stripe_minus_up <- sf::st_sf(this_stripe_minus_up)
        # if there are any other isolines crossing
        all_clines_sf_collapsed <- do.call(rbind, lapply(all_clines_sf, sf::st_sf))
        whether_cross <- c(sf::st_crosses(all_clines_sf_collapsed, sf::st_union(this_stripe_minus_up), sparse = FALSE))
        if (any(whether_cross)) {
          sub_region_ind <- sf::st_intersects(this_stripe_minus_up, all_clines_sf_collapsed[whether_cross, ], sparse = FALSE)
          sub_region_ind <- which(rowSums(sub_region_ind) > 0L)
          subregions <- lapply(sub_region_ind, function(subr) {
            lwgeom::st_split(this_stripe_minus_up[sub_region_ind, ], 
                             sf::st_combine(all_clines_sf_collapsed[whether_cross, ])) |>
              sf::st_collection_extract("POLYGON")
          })
          this_stripe_minus_up <- this_stripe_minus_up[-sub_region_ind, ]
          this_stripe_minus_up <- rbind(do.call(rbind, subregions), this_stripe_minus_up)
        }
        # check if any ups
        this_stripe_minus_up_code <- sapply(1:nrow(this_stripe_minus_up), function(xx) {
          xx <- this_stripe_minus_up[xx, ]
          inds <- sf::st_intersects(xx, grids_pts_sf)
          avglevel <- sapply(inds, function(ii) mean(grids_pts_sf$density_coi_average[ii]))
          avglevel <- findInterval(avglevel, levs, rightmost.closed = TRUE)
          return(avglevel)
        })
        this_stripe_still_up <- this_stripe_minus_up[this_stripe_minus_up_code >= lev_code, ]
        return(this_stripe_still_up)
      })
      any_still_up <- sapply(stripes_up, nrow)
      if (any(any_still_up > 0L)) {
        stripes_up <- do.call(rbind, stripes_up)
        stripes_up <- sf::st_difference(stripes_up, sf::st_union(areas_down))
      } else { stripes_up <- NULL }
    } else { stripes_up <- NULL }
    
    if (!is.null(stripes_up)) {
      areas <- sf::st_union(areas, stripes_up)
    }
  }
  
  # check if there are any other up regions
  if (any(cline_touches_bound)) {
    canvas_minus_areas <- sf::st_difference(canvas_sf, sf::st_union(areas))
    if (!is.null(areas_down)) canvas_minus_areas <- sf::st_difference(canvas_minus_areas, sf::st_union(areas_down))
    if (nrow(canvas_minus_areas) > 0L) {
      # if there are any other isolines crossing
      all_clines_sf_collapsed <- do.call(rbind, lapply(all_clines_sf, sf::st_sf))
      whether_cross <- c(sf::st_crosses(all_clines_sf_collapsed, sf::st_union(canvas_minus_areas), sparse = FALSE))
      if (any(whether_cross)) {
        sub_region_ind <- sf::st_intersects(canvas_minus_areas, all_clines_sf_collapsed[whether_cross, ], sparse = FALSE)
        sub_region_ind <- which(rowSums(sub_region_ind) > 0L)
        subregions <- lapply(sub_region_ind, function(subr) {
          lwgeom::st_split(canvas_minus_areas[sub_region_ind, ], 
                           sf::st_combine(all_clines_sf_collapsed[whether_cross, ])) |>
            sf::st_collection_extract("POLYGON")
        })
        canvas_minus_areas <- canvas_minus_areas[-sub_region_ind, ]
        canvas_minus_areas <- rbind(do.call(rbind, subregions), canvas_minus_areas)
      }
      # check if any ups
      canvas_minus_areas_code <- sapply(1:nrow(canvas_minus_areas), function(xx) {
        xx <- canvas_minus_areas[xx, ]
        inds <- sf::st_intersects(xx, grids_pts_sf)
        if (length(unlist(inds)) == 0L) {
          inds <- sf::st_intersects(sf::st_buffer(xx, dist = mean(spe@metadata$grid_info$xstep, spe@metadata$grid_info$ystep)/2), grids_pts_sf)
        }
        avglevel <- sapply(inds, function(ii) mean(grids_pts_sf$density_coi_average[ii]))
        avglevel <- findInterval(avglevel, levs, rightmost.closed = TRUE)
        return(avglevel)
      })
      if (any(canvas_minus_areas_code >= lev_code)) {
        canvas_minus_areas_still_up <- canvas_minus_areas[canvas_minus_areas_code >= lev_code, ]
        areas <- sf::st_union(areas, canvas_minus_areas_still_up)
      }
    }
  }
  
  areas <- sf::st_as_sf(sf::st_union(areas))
}

