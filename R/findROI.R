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
#' @param zoom.in Logical. For very large ROIs, whether to zoom in and try
#' to get more refined ROIs. 
#' @param zoom.in.size A numeric scaler. Smallest size of an ROI to be able
#' to zoom in. Default is 500L. 
#' @param ... Other parameters that passed to walktrap.community when method =
#' "walktrap".
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
                    zoom.in = FALSE, zoom.in.size = 500L,
                    ...) {

  grid_data <- spe@metadata$grid_density
  grid_type <- spe@metadata$grid_info$grid_type

  coi_clean <- `if`(is.null(coi),"overall",cleanName(coi))
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

    eps <- `if` (grid_type=="hex",
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
    g <- make_graph_new(grid_data_filter$node_x,
                     grid_data_filter$node_y,
                     graph_type = {
                       if(grid_type=="hex") "hex"
                       else if (diag.nodes) "diag"
                       else "square"
                     })

    w <- grid_data_filter$density_coi_average
    igraph::E(g)$weight <- (w[igraph::tail_of(g,igraph::E(g))] +
                         w[igraph::head_of(g,igraph::E(g))])/2
    
    if (diag.nodes) {
      is_diag <- !is.na(igraph::E(g)$diag)
      igraph::E(g)$weight[is_diag] <- igraph::E(g)$weight[is_diag]/sqrt(2)
    }
    igraph::V(g)$name <- grid_data_filter$node
    
    g_community <- switch (method,
                           walktrap = igraph::cluster_walktrap(g, ...),
                           connected = igraph::groups(igraph::components(g)),
                           eigen = igraph::cluster_leading_eigen(g),
                           greedy = igraph::cluster_fast_greedy(g))
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
    seq_along(g_community),
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
    component_list$xcoord <- sort(unique00(spe@metadata$grid_density$x_grid))[
      as.numeric(component_list$x)*2 - as.numeric(component_list$y)%%2
    ]
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
                                      labels = seq_along(unique(rois_filtered$component))
    )
  }
  
  coi_clean_output <- paste(c(coi_clean,"roi"), collapse="_")
  # spe@metadata$coi <- coi
  # spe@metadata$ngrid.min <- ngrid.min
  rois_filtered <- S4Vectors::DataFrame(rois_filtered)
  S4Vectors::metadata(rois_filtered) <- list(densities = dens_cols)
  spe@metadata[[coi_clean_output]] <- rois_filtered

  return(spe)
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
make_graph_new <- function(node_x,node_y,graph_type = c("square","diag","hex")) {
  # hex & igraph both go row-by-row but spatstat goes col-by-col
  if(graph_type!="hex") {
    temp<-node_x
    node_x<-node_y
    node_y<-temp
  }
  range_x <- range(node_x)
  range_y <- range(node_y)
  nx <- diff(range_x)+1
  ny <- diff(range_y)+1
  g <- igraph::make_lattice(dimvector=c(nx,ny))
  # Add edges as needed
  switch(match.arg(graph_type),
         hex = {
           new_edges <- rep(seq_len(nx*(ny-1))[-(nx*seq_len(ny-1))],each = 2)
           if(range_y[1]%%2==1){
             new_edges <- new_edges + rep_len(c(rep.int(c(1,nx),nx-1),
                                                rep.int(c(0,nx+1),nx-1)),
                                              length(new_edges))
           } else{
             new_edges <- new_edges + rep_len(c(rep.int(c(0,nx+1),nx-1),
                                                rep.int(c(1,nx),nx-1)),
                                              length(new_edges))
           }
           g <- igraph::add_edges(g, new_edges)
         },
         diag = {
           new_edges <- rep(seq_len(nx*(ny-1))[-(nx*seq_len(ny-1))],each = 4)
           new_edges <- new_edges + c(0, nx+1, 1, nx)
           g <- igraph::add_edges(g, new_edges, attr = list(diag = TRUE))
         })
  ## Remove filtered vertices
  # This change the row/col into index compatible with vertices in g
  keep_vertices <- node_x - range_x[1] + 1 + (node_y - range_y[1])*nx
  g <- g - seq_len(length(g))[-keep_vertices]
  return(g)
}