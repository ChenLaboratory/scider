#' Build a niche assay based on the profile of neighbouring cells
#'
#' @param spe A SpatialExperiment object
#' @param at Option of cell or grid neighbourhood
#' @param nbrs_name Name of the neighbour list in spe@metadata$grid[[at]]
#' @param group.by Character vector to group neighbours cell by. Should be in 
#' either colData(spe) or spe@metadata$grid_density, depending 
#' on "at". Multiple groups can be used. See details
#' @param use_weight Whether to scale each nbr based on its weight
#' @details
#' For numerical group, result will be sum of nbrs for each cell. For
#' categorical group (factor/string), result will be counts of nbrs belonging in 
#' category
#' @return A matrix where rows are cells/grid points and cols are groups based 
#' on group.by 
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe <- findNbrsSpatial(spe,k=30)
#' niche = getNiche(spe,at="cell",group.by="cell_type")
getNiche <- function(spe,
                     at=c("cell","grid"),
                     nbrs_name = NULL,
                     group.by,
                     use_weight = FALSE) {
  at = match.arg(at)
  # some check for group.by + meaningful error messages
  if (at=="cell") {
    groups <- spe@colData[group.by]
    # row_names <- colnames(spe)
  } else { # grid
    groups <- spe@metadata$grid_density[group.by]
    # row_names <- rownames(spe@metadata$grid_density)
  }
  groups = as.data.frame(groups)
  n_row = nrow(groups)
  groups = cbind(seq_len(n_row),groups)
  
  # Convert categorical/factor/string groups to wide matrix
  groups_matrix = matrix(ncol=0,nrow=n_row,dimnames=list(rownames(groups)))
  for (i in 2:ncol(groups)) {
    if (is.numeric(groups[[i]])) {
      groups_matrix = cbind(groups_matrix,groups[[i]])
    } else {
      groups_matrix = cbind(groups_matrix,+(table(groups[c(1,i)])!=0))
    }
  }
  
  nbrs <- spe@metadata$nbrs[[at]][[nbrs_name%||%length(spe@metadata$nbrs[[at]])]]
  # Aggregate each group for cells based on their nbrs
  for (g in 1:ncol(groups_matrix)) {
    group = groups_matrix[,g]
    if (use_weight) {
      niche <- lapply(seq_along(nbrs$index),function(i){
        sum(group[nbrs$index[[i]]]*nbrs$weight[[i]])
      })
    } else {
      niche <- lapply(nbrs$index,function(i){
        sum(group[i])
      })
    }
    groups_matrix[,g] = unlist(niche)
  }
  
  return(groups_matrix)
}
