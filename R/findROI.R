#' Find ROIs based on cell type-specific densities via graph-based method.
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' Default to all cell types.
#' @param probs A numeric scalar. The threshold of proportion that used to
#'  filter grids by density. Default to 0.85.
#' @param min.density A numeric value. The cut-off value used to filter grids
#' by density. Default is NULL and overwrites probs.
#' @param ngrid.min An integer. The minimum number of grids required for
#' defining a ROI. Default to 20.
#' @param method The community dectection method to be used, possible options
#' are greedy, walktrap, connected, hdbscan, eigen or dbscan. 
#' Default to greedy, can be abbreviated.
#' @param diag.nodes Logical. Set this to TRUE to allow diagonal grid points
#' to be adjacent nodes.
#' @param sequential.roi.name Logical. Set this to FALSE if you want the
#' original ROI name before
#' filtering are retained.
#' @param directed Logical. For graph-based approaches, whether to build a
#' directed graph. 
#' @param zoom.in Logical. For very large ROIs, whether to zoom in and try
#' to get more refined ROIs. 
#' @param zoom.in.size A numeric scaler. Smallest size of an ROI to be able
#' to zoom in. Default is 500L. 
#' @param ... Other parameters that passed to walktrap.community when method =
#' "walktrap".
#' @param fast Whether to create graph using the newer method. May give 
#' different clustering to stochastic methods like walktrap and greedy 
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @importFrom dbscan hdbscan
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
findROI <- function(spe, coi = NULL,
                    probs = 0.85,
                    min.density = NULL, 
                    ngrid.min = 20,
                    method = c("greedy", "walktrap", "connected", "hdbscan", "eigen", "dbscan"),
                    diag.nodes = FALSE,
                    sequential.roi.name = TRUE, 
                    directed = FALSE, 
                    zoom.in = FALSE, zoom.in.size = 500L, 
                    fast = FALSE, ...) {

  grid_data <- spe@metadata$grid_density
  grid_type <- spe@metadata$grid_info$grid_type

  if (is.null(coi)) coi <- "overall"
  if (length(coi) >= 2) coi <- coi[coi!="overall"]
  coi_clean <- janitor::make_clean_names(coi)
  dens_cols <- paste("density", coi_clean, sep = "_")
  
  if (!all(dens_cols %in% colnames(grid_data))) {
    stop("Density of COI is not yet computed.")
  }
  
  method <- match.arg(method)

  grid_data$density_coi_average <- rowSums(as.matrix(grid_data[, which(colnames(grid_data) %in% dens_cols), drop = FALSE]))
  
  if (!is.null(min.density)) {
    message("Overwriting the probs argument. Grids are filtered by the min.density value. ")
    kp <- grid_data$density_coi_average >= min.density
  } else {
    kp <- grid_data$density_coi_average >= quantile(grid_data$density_coi_average, probs = probs)
  }
  grid_data_filter <- grid_data[kp, ]
  
  # clustering approach
  if (identical(method, "hdbscan")) {
    message(paste("For hdbscan, using minPts = ", ngrid.min, sep = ""))
    cl <- dbscan::hdbscan(grid_data_filter[, c("x_grid", "y_grid")], minPts = ngrid.min) #, ...)
    
    if(max(cl$cluster)==0) {
      stop("No clusters detected. Try a different method.")
    }
    
    cls <- setdiff(sort(unique(cl$cluster)), 0)
    g_community <- lapply(cls, function(cc) { grid_data_filter$node[cl$cluster == cc] })
  } else if (identical(method, "dbscan")) {
    message(paste("For dbscan, using minPts = ", ngrid.min, sep = ""))

    eps = `if` (grid_type=="hex",
                2*diff(spe@metadata$grid_info$xlim)/spe@metadata$grid_info$xbins,
                sum(spe@metadata$grid_info$xstep, spe@metadata$grid_info$ystep))
    cl <- dbscan::dbscan(grid_data_filter[, c("x_grid", "y_grid")], 
                         eps = eps, 
                         weights = grid_data_filter[["density_coi_average"]],
                         minPts = ngrid.min)
    
    if(max(cl$cluster)==0) {
      stop("No clusters detected. Try a different method.")
    }
    
    cls <- setdiff(sort(unique(cl$cluster)), 0)
    g_community <- lapply(cls, function(cc) { grid_data_filter$node[cl$cluster == cc] })
  } else {
    # network approaches
    if (!fast) {
    if (grid_type=="hex") {
      adj_edges <- do.call(rbind, lapply(seq_len(nrow(grid_data_filter)), function(ii) {
        adjacent_grids_hex(spe, grid_data_filter$node_x[ii], grid_data_filter$node_y[ii])}))
    } else {
      if (diag.nodes) {
        adj_edges <- do.call(rbind, lapply(seq_len(nrow(grid_data_filter)), function(ii) {
          adjacent_grids_with_corner(spe, grid_data_filter$node_x[ii], grid_data_filter$node_y[ii])}))
        } else {
          adj_edges <- do.call(rbind, lapply(seq_len(nrow(grid_data_filter)), function(ii) {
            adjacent_grids(spe, grid_data_filter$node_x[ii], grid_data_filter$node_y[ii])}))
        }
      }
    keep <- (adj_edges$node1 %in% grid_data_filter$node) & (adj_edges$node2 %in% grid_data_filter$node)
    adj_edges <- adj_edges[keep, ]
    # not allowing self-connected nodes
    keep <- adj_edges$node1 != adj_edges$node2
    adj_edges <- adj_edges[keep, ]
    adj_edges$node1_wt <- grid_data_filter$density_coi_average[match(adj_edges$node1, grid_data_filter$node)]
    adj_edges$node2_wt <- grid_data_filter$density_coi_average[match(adj_edges$node2, grid_data_filter$node)]
    adj_edges$weight <- (adj_edges$node1_wt + adj_edges$node2_wt) * 0.5
    
    if (diag.nodes && grid_type=="square") {
      adj_edges$weight[adj_edges$class == "corner"] <- adj_edges$weight[adj_edges$class == "corner"] / sqrt(2)
    }
    
    df_edges <- adj_edges[, c("node1", "node2", "weight")]

    g <- igraph::graph_from_data_frame(df_edges, directed = directed)
    } else {
      g = make_graph_new(grid_data_filter$node_x,
                         grid_data_filter$node_y,
                         graph_type = {
                           if(grid_type=="hex") "hex"
                           else if (diag.nodes) "diag"
                           else "square"
                         })
      
      # sorted order for grid_data_filter to match g
      order = order(grid_data_filter$node_y, grid_data_filter$node_x)
      
      w = grid_data_filter[order, "density_coi_average"]
      igraph::E(g)$weight = (w[igraph::tail_of(g,igraph::E(g))] + 
                             w[igraph::head_of(g,igraph::E(g))])/2
      
      if (diag.nodes) {
        is_diag = !is.na(igraph::E(g)$diag)
        igraph::E(g)$weight[is_diag] = igraph::E(g)$weight[is_diag]/sqrt(2)
      }
      igraph::V(g)$name = grid_data_filter[order, "node"]
    }


    if (method == "walktrap") {
      g_community <- igraph::cluster_walktrap(g, ...)
    }
    if (method == "connected") {
      g_community <- igraph::groups(igraph::components(g))
    }
    if (method == "eigen") {
      g_community <- igraph::cluster_leading_eigen(g)
    }
    if (method == "greedy") {
      if (!fast) {
        df_edges <- adj_edges[adj_edges$node1_wt < adj_edges$node2_wt, c("node1", "node2", "weight")]
        g <- igraph::graph_from_data_frame(df_edges, directed = directed)
      }
      g_community <- igraph::cluster_fast_greedy(g)
    }
    if (zoom.in) {
      connected_groups <- igraph::groups(igraph::components(g))
      g_community <- lapply(names(connected_groups), function(ind) {
        this_grp <- connected_groups[[ind]]
        if (length(this_grp) > zoom.in.size) {
          subg <- igraph::induced_subgraph(g, this_grp)
          suppressWarnings(subg_community <- igraph::cluster_leading_eigen(subg))
          subg_community <- igraph::communities(subg_community)
          names(subg_community) <- paste(ind, names(subg_community), sep = "-")
        } else {
          subg_community <- list(this_grp)
          names(subg_community) <- ind
        }
        return(subg_community)
      })
      g_community <- do.call(c, g_community)
    }
    
  }
  
  component_list <- do.call(rbind, lapply(
    seq_len(length(g_community)),
    function(ii) {
      data.frame(component = ii, members = g_community[[ii]])
    }
  ))
  
  component_list <- cbind(
    component_list,
    do.call(rbind, strsplit(component_list$members, split = "-"))
  )
  colnames(component_list)[3:4] <- c("x", "y")
  if (grid_type=="hex") {
    temp=sort(unique00(spe@metadata$grid_density$x_grid))
    temp=temp[seq.int(1,length(temp),2)]
    component_list$xcoord <- temp[
      as.numeric(component_list$x)
    ] + (diff(temp[1:2])/2)*(!as.numeric(component_list$y)%%2)
    component_list$ycoord <- sort(unique00(spe@metadata$grid_density$y_grid))[
      as.numeric(component_list$y)
    ]
  } else {
    component_list$xcoord <- spe@metadata$grid_info$xcol[
      as.numeric(component_list$x)
    ]
    component_list$ycoord <- spe@metadata$grid_info$yrow[
      as.numeric(component_list$y)
    ]
  }
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

# Get all neighbouring hexs, even those with negative or out-of-bound index
adjacent_grids_hex = function(spe, node_xx, node_yy) {
  center_cell <- paste(node_xx, node_yy, sep = "-")
  is_odd = node_yy%%2
  cells <- paste(
    c(node_xx-1,node_xx+1,node_xx-is_odd,node_xx+1-is_odd,node_xx-is_odd,node_xx+1-is_odd),
    c(node_yy,node_yy,node_yy+1,node_yy+1,node_yy-1,node_yy-1),
    sep="-"
  )
  return(
    data.frame(
      node1 = rep(center_cell,6),
      node2 = cells
    )
  )
}

# Faster graph maker
#
# @param node_x integer vector of node column
# @param node_y integer vector of node row. Must be same length as node_x
# @param grid_type options of 'square' for square lattice, 'diag' for square lattice 
# with diagonal, and 'hex' for hexagonal lattice.
# 
# Names of vertices are sequential, going row by row:
#     4| 13  14  15  16
#     3| 9   10  11  12
# row 2| 5   6   7   8
#     1| 1   2   3   4
#      ---------------
#        1   2   3   4
#             col 
make_graph_new = function(node_x,node_y,graph_type = c("square","diag","hex")) {
  range_x = range(node_x)
  range_y = range(node_y)
  nx = diff(range_x)+1
  ny = diff(range_y)+1
  g = igraph::make_lattice(dimvector=c(nx,ny))
  
  switch(match.arg(graph_type),
         hex = {
           new_edges = rep(seq_len(nx*(ny-1))[-(nx*seq_len(ny-1))],each = 2)
           if(range_y[1]%%2==1){
             new_edges = new_edges + rep_len(c(rep.int(c(1,nx),nx-1),
                                               rep.int(c(0,nx+1),nx-1)),
                                             length(new_edges))
           } else{
             new_edges = new_edges + rep_len(c(rep.int(c(0,nx+1),nx-1),
                                               rep.int(c(1,nx),nx-1)),
                                             length(new_edges))
           }
           g = igraph::add_edges(g, new_edges)
         },
         diag = {
           new_edges = rep(seq_len(nx*(ny-1))[-(nx*seq_len(ny-1))],each = 4)
           new_edges = new_edges + c(0, nx+1, 1, nx)
           g = igraph::add_edges(g, new_edges, attr = list(diag = TRUE))
         })
  
  ## Remove filtered vertices 
  # This change the row/col into index compatible with vertices in g
  keep_vertices = node_x - range_x[1] + 1 + (
    node_y - range_y[1])*g$dimvector[1]
  g = g - seq_len(length(g))[-keep_vertices]
  return(g)
}