#' Plot cell type composition in each density level of cell of interest.
#'
#' @param spe A SpatialExperiment object.
#' @param contour A character vector of cell type(s) on which the 
#' contour density level is calculated. If NULL, it looks for 
#' 'overall_contour' in colData(spe). Default to NULL.
#' @param id A character. The name of the column of colData(spe)
#' containing the cell type identifiers. Set to 'cell_type' by default.
#' @param roi Character. The name of the group or cell type on which
#' the roi is computed. Default is NULL for no plotting by ROI
#' @param self.included Logical. Whether to include all the cell types in the plot. 
#' Default to TRUE. If FALSE, the cell types specified in 'contour' will not
#' be included in the plot.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#' coi <- "Breast cancer"
#' spe <- gridDensity(spe, coi = coi)
#' spe <- findROI(spe, coi = coi)
#' spe <- getContour(spe, coi = coi)
#' spe <- allocateCells(spe, contour = coi)
#' plotCellCompo(spe, contour = "Breast cancer")
#' plotCellCompo(spe, contour = "Breast cancer", roi = coi)
#'
plotCellCompo <- function(spe, 
                          contour = NULL, 
                          id = "cell_type",
                          roi = NULL, 
                          self.included = TRUE) {
    dat <- as.data.frame(colData(spe))
    
    #    if (!(id %in% colnames(dat)) | !(level.name %in% colnames(dat))) {
    #        stop("Input 'id' and 'level.name' are expected to be in
    #        the column names of the colData of spe.")
    #    }
    
    if (!(id %in% colnames(dat))) {
        stop("Input 'id' is expected to be in the column names of the colData(spe).")
    }
    
    if (is.null(contour)) contour <- "overall"
    contour_clean <- janitor::make_clean_names(contour)
    if(length(contour_clean)>1)
        contour_clean <- paste(sort(contour_clean), collapse="_")
    
    #if(! paste0(contour_clean, "_contour") %in% names(spe@metadata)){
    if(! paste0(contour_clean, "_contour") %in% names(colData(spe)))
        stop("Specified contour not detected in colData(spe). Please run 'getContour' and 'allocateCells' first.")
    level.name <- paste0(contour_clean, "_contour")
    
    if (is.null(roi)) {
        dat <- dat[, c(id, level.name), drop=FALSE]
    } else {
        roi <- gsub("_roi$", "", roi)
        roi <- janitor::make_clean_names(roi)
        roi <- paste(c(sort(roi),"roi"), collapse="_")
        if (is.null(dat[[roi]])) {
            stop(paste(
                roi, " is not found in colData of spe. Please run  
                allocateCells() with 'to.roi=TRUE' first."
            ))
        }
        dat <- dat[, c(id, level.name, roi)]
    }
    
    #if(! (contour == "overall" | all(contour %in% dat[[id]]) ) )
    if ( !all(contour %in% dat[[id]]) && any(contour != "overall")) {
        stop("The cell type(s) on which the contour was computed is not found in colData(spe)[[id]].")
    }
    
    if(!self.included && any(contour != "overall")) dat <- dat[! dat[[id]] %in% contour, ]
    
    if (is.null(roi)) {
        toplot <- calc_proportions(dat, level.name, id)
    } else {
        dat <- dat[dat[[roi]] != "no_roi", ]
        
        grouped_data <- split(dat, as.character(dat[[roi]]))
        
        proportions_by_roi <- lapply(grouped_data, function(group) {
            calc_proportions(group, level.name, id, roi)
        })
        
        toplot <- do.call(rbind, proportions_by_roi)
    }
    col.p <- selectColor(length(unique(toplot[[id]])))

    p <- ggplot(toplot, aes(
        x= .data[[level.name]],
        y = Proportion,
        fill = .data[[id]]
    )) +
        geom_bar(stat = "identity") +
        scale_fill_manual("Cell type", values = col.p) +
        ggtitle(paste0("Cell type composition at each density level of ", 
                       paste(contour, collapse=", "))) +
        theme_classic()
    
    if (is.null(roi)) {
        return(p)
    } else {
        p <- p +
            facet_wrap(~roi)
        
        return(p)
    }
}

calc_proportions <- function(x, level.name, id, roi = NULL) {
    cell_count <- table(x[[level.name]], x[[id]])
    cell_prop <- data.frame(prop.table(cell_count, margin = 1))
    colnames(cell_prop) <- c(level.name, id, "Proportion")
    cell_prop <- cell_prop[!is.na(cell_prop$Proportion), ]
    if (!is.null(roi)) {
        cell_prop$roi <- unique(x[[roi]])
    }
    return(cell_prop)
}

utils::globalVariables(c("Proportion", "cell_type"))
