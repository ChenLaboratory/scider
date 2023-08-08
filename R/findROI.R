#' Find ROIs based on cell type-specific densities via graph-based method.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs). 
#' @param probs A numeric scalar. The threshold of proportion that used to
#'  filter grid by density. Default to 0.85.
#' @param method The community dectection method to be used, either walktrap or connected. Default to walktrap.
#' @param ... Other parameters that passed to walktrap.community.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#' 
#' data("xenium_bc_spe")
#' 
#' coi <- c("Breast cancer", "Fibroblasts")
#' 
#' spe <- gridDensity(spe, coi = coi)
#' 
#' spe <- findROI(spe, coi = coi, method = "walktrap")
#' 

findROI <- function(spe, coi, 
                    probs = 0.85, method = "walktrap", ...) {

  grid_data <- spe@metadata$grid_density

  coi_clean <- janitor::make_clean_names(coi)
  dens_cols <- paste("density", coi_clean, sep="_")
    
  if(!all(dens_cols %in% colnames(grid_data)))
    stop("Density of COI is not yet computed.")

  if (!(method %in% c("walktrap", "connected"))){
    stop("The method chosen is not supported, please choose walktrap or connected.")
  }

  grid_data$density_coi_average <- rowMeans(grid_data[, which(colnames(grid_data) %in% dens_cols), drop=FALSE])

  # Filter grids for COI(s)
  #  kp <- rep(FALSE, nrow(grid_data))
  #  for(i in 1:length(dens_cols)){
  #    sel <- grep(dens_cols[i], colnames(grid_data))
  #   kp <- kp | grid_data[, sel] >= quantile(grid_data[, sel], probs=probs[i])
  #  }

  # Filter grids for COI(s) using average density
  kp <- grid_data$density_coi_average >= quantile(grid_data$density_coi_average, probs=probs)

  grid_data_filter <- grid_data[kp, ]
  
  adj_edges <- do.call(rbind, lapply(1:nrow(grid_data_filter), function(ii) {
    adjacent_grids(grid_data_filter$node_x[ii], grid_data_filter$node_y[ii])
  }))

  keep <- (adj_edges$node1 %in% grid_data_filter$node) & (adj_edges$node2 %in% grid_data_filter$node)
  adj_edges <- adj_edges[keep, ]
  #edge_wts <- (adj_nodes$node_levels[match(adj_edges$node1, adj_nodes$node)] + adj_nodes$node_levels[match(adj_edges$node2, adj_nodes$node)])/2
  edge_wts <- (grid_data_filter$density_coi_average[match(adj_edges$node1, grid_data_filter$node)] + 
               grid_data_filter$density_coi_average[match(adj_edges$node2, grid_data_filter$node)])/2
  
  # create graph
  g <- igraph::graph_from_data_frame(adj_edges, directed = FALSE)
  if (method == "walktrap") {
    g_community <- igraph::cluster_walktrap(g, weights = edge_wts, ...)
  }
  if (method == "connected") {
    g_community <- igraph::groups(igraph::components(g))
  }
  
  #component_list <- lapply(1:length(g_community), function(ii) {
  #  data.frame(component = ii, members = g_community[[ii]])
  #}) %>% do.call(rbind, .)
  
  component_list <- do.call(rbind, lapply(1:length(g_community), function(ii) {
    data.frame(component = ii, members = g_community[[ii]])
  }))
  
  component_list <- cbind(component_list, 
                          do.call(rbind, strsplit(component_list$members, split = "-")))
  colnames(component_list)[3:4] <- c("x", "y")
  component_list$xcoord <- spe@metadata$grid_info$xcol[as.numeric(component_list$x)]
  component_list$ycoord <- spe@metadata$grid_info$yrow[as.numeric(component_list$y)]
  component_list$component <- as.factor(component_list$component)
  
  spe@metadata$components <- component_list
  spe@metadata$coi <- coi
  return(spe)
}


# get edge list
adjacent_grids <- function(node_xx, node_yy) {
  center_cell <- paste(node_xx, node_yy, sep = "-")
  cells <- paste(c(node_xx, node_xx, pmax(node_xx - 1, 1), 
                   pmin(node_xx + 1, spe@metadata$grid_info$dims[1])), 
                 c(pmax(node_yy - 1, 1), 
                   pmin(node_yy + 1, spe@metadata$grid_info$dims[2]), node_yy, node_yy), 
                 sep = "-")
  cells <- cells[cells != center_cell]
  data.frame(node1 = rep(center_cell, length(cells) + 1), node2 = c(center_cell, cells))
}