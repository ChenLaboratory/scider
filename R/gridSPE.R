#' Summarize a SpatialExperiment object at grid-level
#'
#' @param spe A SpatialExperiment object.
#' @param cell.count Logical. Whether to obtain the number of cells within 
#' each group identified by the 'id' column in colData(spe). Default to FALSE.
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default.
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
gridSPE <- function(spe, cell.count = FALSE, id = 'cell_type') {
    if (!("grid_density" %in% names(spe@metadata))) {
        stop("Please run gridDensity before using this function.")
    }

    grid_data <- spe@metadata$grid_density[,(1:2)]

    grids = sf::st_sfc(grid2sf(spe))
    # Assign cells to grids
    xy_allcells <- sf::st_as_sf(
      as.data.frame(SpatialExperiment::spatialCoords(spe)), 
      coords = c("x_centroid", "y_centroid")
    )
    message("Assigning cells to grids.")
    overlap_ind <- sf::st_intersects(xy_allcells, grids, sparse = FALSE)

    # Obtain gene counts at the grid level
    assays <- list()
    assays$counts <- as.matrix(SummarizedExperiment::assay(spe,"counts") %*% overlap_ind)

    # Obtain cell type counts at the grid level
    if(cell.count){
        cell_type_names <- names(table(spe@colData[[id]]))
        cell_counts <- matrix(0, length(cell_type_names), ncol(spe))
        rownames(cell_counts) <- cell_type_names
        for(i in 1:length(cell_type_names)){
            cell_counts[i,] <- spe@colData[[id]] == cell_type_names[i]
        }
        cell_counts <- t(as.matrix(cell_counts %*% overlap_ind))
        cell_counts <- cbind(cell_counts, overall=rowSums(cell_counts))
        cell_counts_names <- colnames(cell_counts) <- paste0("cell_counts_", 
            janitor::make_clean_names(colnames(cell_counts)))
        grid_data[, cell_counts_names] <- cell_counts
    }

    spe_out <- SpatialExperiment::SpatialExperiment(assays = assays, 
                                                    colData = grid_data, 
                                                    rowData = SummarizedExperiment::rowData(spe),
                                                    spatialCoordsNames = c("x_grid", "y_grid"))
    
    spe_out@metadata = spe@metadata
    spe_out@metadata$grid_info$gridLevelAnalysis = TRUE

    return(spe_out)
}
