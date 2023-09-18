#' Find ROIs based on cell type-specific densities via graph-based method.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' @param probs A numeric scalar. The threshold of proportion that used to
#'  filter grid by density. Default to 0.85.
#' @param ngrid.min An integer. The minimum number of grids required for
#' defining a ROI. Default to 20.
#' @param method The community dectection method to be used, either walktrap
#' or connected. Default to walktrap.
#' @param ... Other parameters that passed to walktrap.community.
#' @param diag.nodes Logical. Set this to TRUE to allow diagonal grid points
#' to be adjacent nodes.
#' @param sequential.roi.name Logical. Set this to FALSE if you want the
#' original ROI name before
#' filtering are retained.
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
                    probs = 0.85,
                    ngrid.min = 20,
                    method = "walktrap",
                    diag.nodes = FALSE,
                    sequential.roi.name = TRUE, ...) {
  grid_data <- spe@metadata$grid_density

  coi_clean <- janitor::make_clean_names(coi)
  dens_cols <- paste("density", coi_clean, sep = "_")

  if (!all(dens_cols %in% colnames(grid_data))) {
    stop("Density of COI is not yet computed.")
  }

  if (!(method %in% c("walktrap", "connected"))) {
    stop("The method chosen is not supported, please choose walktrap or connected.")
  }

  grid_data$density_coi_average <- rowMeans(
    as.matrix(grid_data[, which(colnames(grid_data) %in% dens_cols),
      drop = FALSE
    ])
  )

  # Filter grids for COI(s) using average density
  kp <- grid_data$density_coi_average >= quantile(grid_data$density_coi_average,
    probs = probs
  )

  grid_data_filter <- grid_data[kp, ]

  if (!diag.nodes) {
    adj_edges <- do.call(rbind, lapply(seq_len(nrow(grid_data_filter)), function(ii) {
      adjacent_grids(spe, grid_data_filter$node_x[ii], grid_data_filter$node_y[ii])
    }))
  } else {
    adj_edges <- do.call(rbind, lapply(seq_len(nrow(grid_data_filter)), function(ii) {
      adjacent_grids_with_corner(
        spe, grid_data_filter$node_x[ii],
        grid_data_filter$node_y[ii]
      )
    }))
  }

  keep <- (adj_edges$node1 %in% grid_data_filter$node) & (adj_edges$node2 %in%
    grid_data_filter$node)
  adj_edges <- adj_edges[keep, ]
  edge_wts <- (grid_data_filter$density_coi_average[match(
    adj_edges$node1,
    grid_data_filter$node
  )] +
    grid_data_filter$density_coi_average[match(
      adj_edges$node2,
      grid_data_filter$node
    )]) / 2
  # if consider corner nodes, down weight by sqrt(2)
  if (diag.nodes) {
    edge_wts[adj_edges$class == "corner"] <- edge_wts[adj_edges$class == "corner"] / sqrt(2)
  }

  # create graph
  g <- igraph::graph_from_data_frame(adj_edges, directed = FALSE)
  if (method == "walktrap") {
    g_community <- igraph::cluster_walktrap(g, weights = edge_wts, ...)
  }
  if (method == "connected") {
    g_community <- igraph::groups(igraph::components(g))
  }

  component_list <- do.call(rbind, lapply(seq_len(length(g_community)), function(ii) {
    data.frame(component = ii, members = g_community[[ii]])
  }))

  component_list <- cbind(
    component_list,
    do.call(rbind, strsplit(component_list$members, split = "-"))
  )
  colnames(component_list)[3:4] <- c("x", "y")
  component_list$xcoord <- spe@metadata$grid_info$xcol[as.numeric(component_list$x)]
  component_list$ycoord <- spe@metadata$grid_info$yrow[as.numeric(component_list$y)]
  component_list$component <- as.factor(component_list$component)

  # filtering ROIs based on ngrid.min
  filtered <- names(which(table(component_list$component) >= ngrid.min))
  rois_filtered <- component_list[component_list$component %in% filtered, ]

  if (sequential.roi.name) {
    rois_filtered$component <- factor(rank(rois_filtered$component),
      labels = seq(length(unique(rois_filtered$component)))
    )
  }

  spe@metadata$coi <- coi
  spe@metadata$ngrid.min <- ngrid.min
  spe@metadata$roi <- S4Vectors::DataFrame(rois_filtered)

  return(spe)
}


# get edge list
adjacent_grids <- function(spe, node_xx, node_yy) {
  center_cell <- paste(node_xx, node_yy, sep = "-")
  cells <- paste(
    c(
      node_xx, node_xx, pmax(node_xx - 1, 1),
      pmin(node_xx + 1, spe@metadata$grid_info$dims[1])
    ),
    c(
      pmax(node_yy - 1, 1),
      pmin(node_yy + 1, spe@metadata$grid_info$dims[2]), node_yy, node_yy
    ),
    sep = "-"
  )
  cells <- cells[cells != center_cell]
  data.frame(
    node1 = rep(center_cell, length(cells) + 1),
    node2 = c(center_cell, cells)
  )
}


adjacent_grids_with_corner <- function(spe, node_xx, node_yy) {
  center_cell <- paste(node_xx, node_yy, sep = "-")
  all_xx <- c(
    pmax(node_xx - 1, 1), node_xx,
    pmin(node_xx + 1, spe@metadata$grid_info$dims[1])
  )
  all_xx <- all_xx[!duplicated(all_xx)]
  all_yy <- c(
    pmax(node_yy - 1, 1), node_yy,
    pmin(node_yy + 1, spe@metadata$grid_info$dims[2])
  )
  all_yy <- all_yy[!duplicated(all_yy)]
  all_cells <- paste(rep(all_xx, each = length(all_yy)),
    rep(all_yy, length(all_xx)),
    sep = "-"
  )
  all_cells <- all_cells[all_cells != center_cell]
  immediate_cells <- paste(
    c(
      node_xx, node_xx, pmax(node_xx - 1, 1),
      pmin(node_xx + 1, spe@metadata$grid_info$dims[1])
    ),
    c(
      pmax(node_yy - 1, 1),
      pmin(node_yy + 1, spe@metadata$grid_info$dims[2]),
      node_yy, node_yy
    ),
    sep = "-"
  )
  immediate_cells <- immediate_cells[immediate_cells != center_cell]
  corner_cells <- setdiff(all_cells, immediate_cells)
  data.frame(
    node1 = rep(center_cell, length(all_cells) + 1),
    node2 = c(center_cell, immediate_cells, corner_cells),
    class = c(
      "center", rep("immediate", length(immediate_cells)),
      rep("corner", length(corner_cells))
    )
  )
}
