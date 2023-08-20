#' Get contour from density
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' @param bins 
#' @param binwidth 
#' @param breaks 
#'
#' @return An sf object of the contour region of the specified level. 
#' @export
#'
#' @examples
getContour <- function(spe, coi, bins = NULL, binwidth = NULL, breaks = NULL) {
  
  if (is.null(spe@metadata$grid_density))
    stop("Have to calculate grid density, run gridDensity() first!")
  
  dens <- spe@metadata$grid_density
  dups <- duplicated(dens[, c("y_grid", "x_grid"), drop = FALSE], fromLast = TRUE)
  dens <- dens[!dups, , drop = FALSE]
  
  coi_clean <- janitor::make_clean_names(coi)
  dens_cols <- paste("density", coi_clean, sep="_")
  
  if(!all(dens_cols %in% colnames(dens)))
    stop("Density of COI is not yet computed.")
  
  if (length(dens_cols) > 1L) {
    message("Plotting contour of average density of input COIs. ")
    dens$density_coi_average <- rowMeans(dens[, which(colnames(dens) %in% dens_cols), drop=FALSE])
  } else {
    dens$density_coi_average <- dens[, dens_cols]
  }
  dens <- dens[, c(1:5, which(colnames(dens) == "density_coi_average"))]
  
  # levels for contour
  if (is.null(bins) && is.null(binwidth) && is.null(breaks)) {
    warning("Using bins = 10 to draw contours!")
    bins <- 10L
  }
  if (!is.null(bins)) binwidth <- breaks <- NULL
  if (is.null(bins) && !is.null(binwidth)) breaks <- NULL
  
  # filter out negative densities when calculating contours
  dens <- dens[dens$density_coi_average > 0L, ]
  
  # note that when calculating contours, density is not filtered at any quantile cutoff!
  spe@metadata$contour <- compute_group(dens,
                                        z.range = range(dens$density_coi_average, na.rm = TRUE, finite = TRUE),
                                        bins = bins,
                                        binwidth = binwidth,
                                        breaks = breaks,
                                        na.rm = FALSE)
  spe@metadata$contour$level_factor <- as.factor(as.numeric(as.factor(spe@metadata$contour$level)))
  names(spe@metadata)[which(names(spe@metadata) == "contour")] <- paste(coi_clean, "contour", sep = "_")
  
  return(spe)
  
}



#### lower level functions for computing the contours. NEED TO CLEAN UP LATER!!!!
# compute contour groups (grabbed from ggplot2)
xyz_to_isolines <- function(data, breaks) {
  isoband::isolines(
    x = sort(ggplot2:::unique0(data$x_grid)),
    y = sort(ggplot2:::unique0(data$y_grid)),
    z = isoband_z_matrix(data),
    levels = breaks
  )
}

isoband_z_matrix <- function(data) {
  # Convert vector of data to raster
  x_pos <- as.integer(factor(data$x_grid, levels = sort(ggplot2:::unique0(data$x_grid))))
  y_pos <- as.integer(factor(data$y_grid, levels = sort(ggplot2:::unique0(data$y_grid))))
  nrow <- max(y_pos)
  ncol <- max(x_pos)
  raster <- matrix(NA_real_, nrow = nrow, ncol = ncol)
  raster[cbind(y_pos, x_pos)] <- data$density_coi_average
  raster
}

iso_to_path <- function(iso, group = 1) {
  lengths <- vapply(iso, function(x) length(x$x), integer(1))
  
  if (all(lengths == 0)) {
    cli::cli_warn("{.fn stat_contour}: Zero contours were generated")
    return(ggplot2:::data_frame0())
  }
  
  levels <- names(iso)
  xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
  ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
  ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
  item_id <- rep(seq_along(iso), lengths)
  
  # Add leading zeros so that groups can be properly sorted
  groups <- paste(group, sprintf("%03d", item_id), sprintf("%03d", ids), sep = "-")
  groups <- factor(groups)
  
  ggplot2:::data_frame0(
    level = rep(levels, lengths),
    x = xs,
    y = ys,
    piece = as.integer(groups),
    group = groups,
    .size = length(xs)
  )
}

compute_group <-  function(data, z.range, bins = NULL, binwidth = NULL, breaks = NULL, na.rm = FALSE) {
  breaks <- ggplot2:::contour_breaks(z.range, bins, binwidth, breaks)
  isolines <- xyz_to_isolines(data, breaks)
  path_df <- iso_to_path(isolines, data$group[1])
  path_df$level <- as.numeric(path_df$level)
  path_df$nlevel <- scales::rescale_max(path_df$level)
  path_df
}
