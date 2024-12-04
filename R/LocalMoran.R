#' Calculate local Moran for 1 to 2 variables.
#'
#' @param spe A SpatialExperiment object.
#' @param data1 Numeric vector 1. Must be same length as nrow(spe@metadata$grid_density).
#' @param data2 Numeric vector 2 for bivariate local Moran. Must be same length as data1.
#' @param diag.nodes Logical. Set this to TRUE to allow diagonal grid points
#' to be adjacent nodes if grid is square.
#' @param significance_cutoff Cutoff for p-value to filter non-significant clusters
#' @param permutations Number of permutations for p-value.
#' @param seed Integer. For random permutations.
#' @return List of lisa_value, clusters, and pseudo p-value.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe <- gridDensity(spe)
#'
#' plotDensity(spe, coi = "Breast cancer")
#'
#' plotDensity(spe, coi = "Fibroblasts")
#'@export
LocalMoran <- function(spe, data1, data2=data1,diag.nodes=F,
                        significance_cutoff = 0.05,permutations = 999,seed = 123456789) {
  stopifnot(length(data1)==length(data2),
            length(data1) == nrow(spe@metadata$grid_density),
            is.numeric(data1),
            is.numeric(data2))
  
  grid_data <- spe@metadata$grid_density
  g=make_graph_new(grid_data$node_x,
                   grid_data$node_y,
                   graph_type = {
                     if(spe@metadata$grid_info$grid_type=="hex") "hex"
                     else if (diag.nodes) "diag"
                     else "square"
                   })
  data1 = data1[order(grid_data$node_y, grid_data$node_x)]
  data2 = data2[order(grid_data$node_y, grid_data$node_x)]

  # TODO: see if can remove this
  n_nbrs = igraph::degree(g)
  
  nbrs = igraph::adjacent_vertices(g,seq_along(g))
  return(.Call("C_localMoran",nbrs,n_nbrs,data1,data2,significance_cutoff,permutations,seed))
}
