#' Summarize a SpatialExperiment object at grid-level
#'
#' @param spe A SpatialExperiment object.
#' @param cell.count Logical. Whether to obtain the number of cells within 
#' each group identified by the 'id' column in colData(spe). Default to FALSE.
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to 'cell_type' by default.
#' @param split.count.by A character. The name of the column of colData(spe). 
#' When it is not NULL, a grid-level count matrix is calculated for each member 
#' specified in that column of colData(spe) and stored in the assays(spe).
#' Set to 'cell_type' by default. 
#' 
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe <- gridDensity(spe)
#'
#' spe_grid <- gridSPE(spe)
#'
gridSPE <- function(spe, cell.count = FALSE, id = 'cell_type', split.count.by = id) {
  if (!("grid_density" %in% names(spe@metadata))) {
    stop("Please run gridDensity before using this function.")
  }
  
  grid_data <- spe@metadata$grid_density[,(1:5)]
  rownames(grid_data) <- grid_data$node
  grid_info <- spe@metadata$grid_info
  assay_matrix <- spe@assays@data[["counts"]]

  #Get vector of which polygon each cell belong to
  xy_allcells <- SpatialExperiment::spatialCoords(spe)
  if(grid_info$grid_type=="hex") {
    poly_index <- hexDensity::xy2hcell(x=xy_allcells[,1],y=xy_allcells[,2],
                                       xbins=grid_info$xbins,
                                       xbnds=grid_info$xlim,
                                       ybnds=grid_info$ylim,
                                       shape=grid_info$shape)
  } else {
    rr <- (xy_allcells[,1]-grid_info$xlim[1])%/%grid_info$xstep
    rr <- pmin.int(rr,grid_info$dim[1]-1)
    cc <- (xy_allcells[,2]-grid_info$ylim[1])%/%grid_info$ystep
    cc <- pmin.int(cc,grid_info$dim[2]-1)
    poly_index <- (rr*grid_info$dims[2]+cc+1)
  }
  
  # Obtain gene counts at the grid level
  assays <- list()
  if(is.null(split.count.by)){
    assays$counts <- matrix(data=0,
                            nrow=nrow(spe),ncol=nrow(spe@metadata$grid_density),
                            dimnames = list(rownames(spe)))
    
    aggCounts <- colsum(assay_matrix, poly_index, reorder = FALSE)
    assays$counts[,unique(poly_index)] <- aggCounts
  } else {
    split_names <- names(table(spe@colData[[split.count.by]]))
    split_names_clean <- janitor::make_clean_names(split_names)
    for(i in seq_along(split_names)){
      sub <- which(spe@colData[[split.count.by]] == split_names[i])
      assays[[i+1]] <- matrix(data=0,
                              nrow=nrow(spe),ncol=nrow(spe@metadata$grid_density),
                              dimnames = list(rownames(spe)))
      aggCounts <- colsum(assay_matrix[,sub], poly_index[sub], reorder = FALSE)
      assays[[i+1]][,unique(poly_index[sub])] <- aggCounts
    }
    assays[[1]] <- Reduce(`+`, assays[-1])
    names(assays) <- c("counts", paste0("counts_", split_names_clean))
  }
  
  # Obtain cell type counts at the grid level
  if(cell.count){
    # Sorting is just for backward consistency
    cell_type_names <- janitor::make_clean_names(sort(unique(spe@colData[[id]])))

    poly_counts <- matrix(0,nrow(spe@metadata$grid_density),length(cell_type_names),
                          dimnames=list(NULL,cell_type_names))
    cell_types <- as.numeric(factor(spe@colData[[id]]))

    for (i in 1:ncol(spe)) {
      poly_counts[poly_index[i],cell_types[i]] <- poly_counts[poly_index[i],cell_types[i]] + 1
    }
    poly_counts <- cbind(poly_counts, overall=rowSums(poly_counts))
    grid_data$cell_count <- S4Vectors::DataFrame(poly_counts)
  }
  grid_data$LibSize <- colSums(assays[[1]])
  grid_data$nGenes <- colSums(assays[[1]]!=0)
  
  spe_out <- SpatialExperiment::SpatialExperiment(assays = assays,
                                                  colData = grid_data,
                                                  rowData = SummarizedExperiment::rowData(spe),
                                                  spatialCoordsNames = c("x_grid", "y_grid"))
  
  spe_out@metadata <- spe@metadata
  spe_out@metadata$grid_info$gridLevelAnalysis <- TRUE
  
  return(spe_out)
}
