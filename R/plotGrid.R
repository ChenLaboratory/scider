#' Plot grid from metadata. 
#'
#' @param spe A SpatialExperiment object.
#' @param group.by values to group polygons by. Must be in colData of spe. 
#' If NULL, will try with col if available.
#' @param cols Colour palette. Can be a vector of colours or a function 
#' that accepts an integer n and return n colours.
#' @param pol.alpha alpha of points between 0 and 1.
#' @param pol.border Boolean. Whether to draw border for each polygon.
#' @param probs Numeric value between 0 and 1, used for filtering
#' uninformative grid. Only applicable for continuous values.
#' 
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe <- gridDensity(spe)
#'
#' plotGrid(spe, group.by = "x_grid")
#'
plotGrid <- function(spe,
                     group.by = NULL,
                     cols = NULL,
                     pol.border = FALSE,
                     pol.alpha = 1,
                     probs=0) {
  grid_data <- as.data.frame(spe@metadata$grid_density)
  
  xstep <- spe@metadata$grid_info$xstep
  ystep <- spe@metadata$grid_info$ystep
  
  # Groups
  if (!is.null(group.by)) {
    group = grid_data[[group.by]]
  } else if (!is.null(cols)) {
    group = factor(rep_len(cols,nrow(grid_data)),levels=unique(cols))
  } else {
    group = NULL
  }
  
  # Filter
  isContinuous = !is.null(group.by) && is.numeric(grid_data[[group.by]])
  if (isContinuous) {
    kp <- group >= quantile(group, probs = probs)
    group <- group[kp]
    grid_data <- grid_data[kp,]
  }
  
  n_colour = length(unique(group))
  if (is.null(cols)&&is.null(group.by)) {
    col.p = NULL
  } else if (is.null(cols)) {
    if (isContinuous) col.p <- col.spec
    else col.p <- selectColor(n_colour)
  } else if (is.function(cols)) {
    col.p <- as.character(cols(n_colour))
  } else {
    col.p <- as.character(cols)
    if (!is.null(group.by) && !isContinuous) col.p <- rep_len(col.p,n_colour)
    else if (!isContinuous) col.p <- rep_len(unique(col.p), n_colour)
  }
  
  
  # Plotting
  p <- ggplot() + 
    geom_sf(
      data = sf::st_as_sfc(grid2sf(spe,
                                   grid_data$node_x,
                                   grid_data$node_y)),
      aes(
        fill = group
      ),alpha=pol.alpha,
      color = if (pol.border) "black" else NA
    ) +
    theme_classic() +
    labs(x = "x", y = "y", fill = group.by) +
    lims(
      x = c(
        min(grid_data[, "x_grid"]) - xstep/2,
        max(grid_data[, "x_grid"]) + xstep/2
      ),
      y = c(
        min(grid_data[, "y_grid"]) - ystep/2,
        max(grid_data[, "y_grid"]) + ystep/2
      )
    )
  
  if (isContinuous) {
    p = p + scale_fill_gradientn(colours = rev(col.p))
  } else {
    p = p + scale_fill_manual(values = col.p)
  }
  return(p)
}
