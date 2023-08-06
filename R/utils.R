# get edge list
adjacent_grids <- function(node_xx, node_yy) {
  center_cell <- paste(node_xx, node_yy, sep = "-")
  cells <- paste(c(node_xx, node_xx, pmax(node_xx - 1, 1), 
                   pmin(node_xx + 1, metadata(spe)$grid_info$dims[1])), 
                 c(pmax(node_yy - 1, 1), 
                   pmin(node_yy + 1, metadata(spe)$grid_info$dims[2]), 
                   node_yy, node_yy), sep = "-")
  cells <- cells[cells != center_cell]
  data.frame(node1 = rep(center_cell, length(cells) + 1), 
             node2 = c(center_cell, cells))
}