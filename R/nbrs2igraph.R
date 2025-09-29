#' Convert nbrs (in spe@metadata$nbrs) into igraph's graph 
#' 
#' 
#' @details
#' nbrs should be a list containing index & weight
#' @export
nbrs2igraph = function(nbrs,directed=FALSE){
  interleaves <- as.vector(
    rbind(rep.int(seq_along(nbrs$index),times=lengths(nbrs$index)),
          unlist(nbrs$index)))
  g <- igraph::make_graph(interleaves,directed=directed) #TODO: check direction
  igraph::E(g)$weight = unlist(nbrs$weight)
  
  g <- igraph::simplify(g,edge.attr.comb = "first")
  return(g)
}