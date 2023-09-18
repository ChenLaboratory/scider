#' Visualising an sf object (for internal use only at the moment)
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of length 1 of the cell type of interest (COIs).
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default.
#' @param overlay Character vector. Either plot overlay on density or cells.
#' @param sub.level Numeric vector of length 1 or 2, identifies which density level
#' to plot. When length is 1, plot the density region above this level. When
#' length is 2, plot the density region between the two levels.
#'
#' @return A ggplot object.
#'
plotContourRegion <- function(spe,
                              coi,
                              id = "cell_type",
                              overlay = c("density", "cell"),
                              sub.level) {
  if (length(coi) > 1L) {
    stop("coi must be of length 1!")
  }

  if (!(coi %in% colData(spe)[[id]])) {
    stop("coi not in colData(spe)[[id]]!")
  }

  if (!(length(sub.level) %in% c(1L, 2L))) {
    stop("sub.level must be of either length 1 or 2!")
  }

  coi_clean <- janitor::make_clean_names(coi)
  coi_clean_contour <- paste(coi_clean, "contour", sep = "_")
  coi_clean_density <- paste("density", coi_clean, sep = "_")

  if (!coi_clean_contour %in% names(spe@metadata)) {
    stop("Contour of coi doesn't exist. Please run getContour() first!")
  }

  contour_data <- as.data.frame(spe@metadata[[coi_clean_contour]])
  levs <- unique(contour_data$level)
  nlevs <- length(levs)

  if (any(!(sub.level %in% levs))) {
    stop("Non-existing level(s) in sub.level! Check unique(spe@metadata$`contour_data`$level).")
  }

  if (length(overlay) == 2) {
    overlay <- "density"
  }

  xlim <- spe@metadata$grid_info$xlim
  ylim <- spe@metadata$grid_info$ylim
  plot.xlim <- xlim + c(-1e-10, 1e-10)
  plot.ylim <- ylim + c(-1e-10, 1e-10)
  w <- spe@metadata$grid_info$xstep
  h <- spe@metadata$grid_info$ystep

  if (length(sub.level) == 1L) {
    density_at_level <- unique(contour_data$cutoff[contour_data$level == sub.level])
  }
  if (length(sub.level) == 2L) {
    density_at_level <- unique(contour_data$cutoff[contour_data$level %in% sub.level])
    density_at_level <- sort(density_at_level, decreasing = FALSE)
  }

  if (overlay == "cell") {
    sub <- grep(coi, colData(spe)[[id]])
    p <- ggplot() +
      geom_point(
        data = as.data.frame(spatialCoords(spe[, sub])),
        aes(x = x_centroid, y = y_centroid), size = 0.02,
        stroke = 0.2, alpha = 0.4
      )
  } else if (overlay == "density") {
    dens_df <- as.data.frame(spe@metadata$grid_density)
    colnames(dens_df)[which(colnames(dens_df) == coi_clean_density)] <- "density"
    if (length(sub.level) == 1L) {
      kp_grids <- dens_df$density >= density_at_level
      kp_isolines <- contour_data$level == sub.level
      lev_name <- sub.level
    }
    if (length(sub.level) == 2L) {
      kp_grids <- dens_df$density >= density_at_level[1] & dens_df$density < density_at_level[2]
      kp_isolines <- contour_data$level %in% sub.level
      lev_name <- paste0(sub.level, collapse = "-")
    }
    p <- ggplot() +
      geom_tile(data = dens_df[kp_grids, ], aes(
        x = x_grid,
        y = y_grid, width = w,
        height = h
      ), fill = "lightgrey") +
      geom_path(data = contour_data[kp_isolines, ], aes(
        x = x, y = y,
        group = group
      )) +
      ggtitle(paste("Level = ", lev_name, sep = ""))
  } else {
    stop("Overlay should either be cell or density.")
  }

  # overlay the sf region
  if (length(sub.level) == 1L) {
    area <- contour2sf(spe,
      contour = coi_clean_contour, coi = coi,
      cutoff = density_at_level
    )
    p <- p + geom_sf(
      data = area, alpha = 0.2, fill = "tomato",
      linewidth = 0.04
    )
  }
  if (length(sub.level) == 2L) {
    area_low <- contour2sf(spe,
      contour = coi_clean_contour,
      coi = coi, cutoff = density_at_level[1]
    )
    area_high <- contour2sf(spe,
      contour = coi_clean_contour,
      coi = coi, cutoff = density_at_level[2]
    )
    area_diff <- sf::st_difference(area_low, area_high)
    if (overlay == "cell") {
      p <- p + geom_sf(
        data = area_diff, alpha = 0.4,
        fill = "tomato", linewidth = 0.2
      )
    }
    if (overlay == "density") {
      p <- p + geom_sf(
        data = area_diff, alpha = 0.2,
        fill = "tomato", linewidth = 0.04
      )
    }
  }

  p <- p +
    scale_x_continuous(limits = plot.xlim) +
    scale_y_continuous(limits = plot.ylim) + labs(x = "x", y = "y")

  return(p)
}

utils::globalVariables(c(
  "x_centroid", "y_centroid", "x_grid",
  "y_grid", "x", "y",
  "group"
))
