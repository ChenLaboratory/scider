#' Plot grid from metadata. 
#'
#' @param spe A SpatialExperiment object.
#' @param reverseY Reverse y coordinates.
#' @param group.by values to group polygons by. Must be in 
#' spe@metadata$grid_density, or colData(spe) if gridLevelAnalysis is TRUE. 
#' If NULL, will try with cols if available.
#' @param feature Feature to group polygons by. Must be in rownames(spe).
#' @param assay Name of assay to use for plotting feature.
#' @param type Transformation to apply for the group/feature. Options are "raw"
#' , "log", "cpm", "logcpm", or a function that accepts and returns a vector of 
#' the same length.
#' @param cols Colour palette. Can be a vector of colours or a function 
#' that accepts an integer n and return n colours.
#' @param pol.alpha alpha of points between 0 and 1.
#' @param pol.border Boolean. Whether to draw border for each polygon.
#' @param probs Numeric value between 0 and 1, used for filtering
#' uninformative grid. Only applicable for continuous values.
#' @param label label for the legend
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
plotGrid <- function(spe, reverseY = FALSE,
                     group.by = NULL,
                     feature = NULL,
                     assay = "counts",
                     type = c("raw","log","cpm","logcpm"),
                     cols = NULL,
                     pol.border = FALSE,
                     pol.alpha = 1,
                     probs = 0,
                     label = NULL) {
  grid_data <- spe@metadata$grid_density
  if(!is.null(spe@metadata$grid_info$gridLevelAnalysis)) {
    grid_data = c(grid_data,colData(spe))
  }
  
  xlim <- range(grid_data[, "x_grid"]) + c(-1,1)*spe@metadata$grid_info$xstep/2
  ylim <- range(grid_data[, "y_grid"]) + c(-1,1)*spe@metadata$grid_info$ystep/2
  
  group = col.p = NULL
  # Groups. Order is: grid_density -> colData -> assays -> cols
  if (!is.null(group.by) && group.by %in% colnames(grid_data)) {
    group = grid_data[[group.by]]
    if (is.null(label)) label = group.by
  } else if (!is.null(feature) && feature %in% rownames(spe)) {
    group = SummarizedExperiment::assays(spe)[[assay]][feature,]
    if (is.null(label)) label = feature
  } else if (!is.null(cols) && !is.function(cols)) {
    group = factor(rep_len(cols,nrow(grid_data)),levels=unique(cols))
    col.p = rep_len(unique(cols), length(unique(cols)))
  }
  
  # Type
  if (is.character(type)) {
    type = switch(match.arg(type),
                  raw = NULL,
                  log = function(x) {log2(x+1)},
                  cpm = function(x) {
                    (x+0.5)/colSums(as.matrix(spe@assays@data[[assay]]))*1e6
                  },
                  logcpm = function(x) {
                    log2((x+0.5)/colSums(as.matrix(spe@assays@data[[assay]]))*1e6)
                  })
  }
  if(is.function(type)) {
    tryCatch({group = type(group)},
             error = function(e){
               message("Error when applying 'type'. Skipping 'type'.")
             })
  }
  isContinuous = is.numeric(group)
  
  # Filter
  if (isContinuous) {
    kp <- group >= quantile(group, probs = probs)
    group <- group[kp]
    grid_data <- grid_data[kp,]
  }
  
  # Colours
  if (!is.null(group) && is.null(col.p)) {
      n_colour = length(unique(group))
      if (is.null(cols)) { # Default palette
          col.p = `if`(isContinuous, col.spec, selectColor(n_colour))
      } else if (is.function(cols)) { # cols is function
          col.p = cols(n_colour)
      } else { # cols is vector
          col.p = `if`(isContinuous, cols, rep_len(cols,n_colour))
      }
  }
  # Plotting
  p <- ggplot() + 
    geom_sf(
      data = sf::st_as_sfc(grid2sf(spe,
                                   grid_data$node_x,
                                   grid_data$node_y,
                                   reverseY=reverseY)),
      aes(
        fill = group
      ),alpha=pol.alpha,
      color = if (pol.border) "black" else NA
    ) +
    theme_classic() +
    labs(x = "x", y = "y", fill = label) +
    lims(
      x = xlim,
      y = ylim
    )
  
  if (isContinuous) {
    p = p + scale_fill_gradientn(colours = rev(col.p))
  } else {
    p = p + scale_fill_manual(values = col.p)
  }
  return(p)
}