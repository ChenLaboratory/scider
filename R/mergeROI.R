#' Manually merge ROIs
#'
#' @param spe A SpatialExperiment object.
#' @param merge.list A (named) list of vectors of ROI ids to be merged.
#' Each vector in the list should be of length greater than or
#' equal to 2. If no name is specified, the merged ROI will be
#' named by concatenating ROIs being merged.
#' @param remove.ids Optional. A vector of ROI ids to be removed. 
#' @param id Character. The name of the column in \code{spe@metadata$roi}
#' that stores the ROIs to be merged. Default is "component".
#' @param rename Logical. If TRUE, names of merge.list are ignored.
#' ROIs will be given a new name. For the unmerged ROIs, their new
#' names are not necessarily the same as those before merging.
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
#' spe <- mergeROI(spe, list("1-2" = 1:2))
#'
mergeROI <- function(spe, 
                     merge.list = NULL, 
                     remove.ids = NULL, 
                     id = "component",
                     rename = FALSE) {
  # check ROI exists
  if (is.null(spe@metadata$roi)) {
    stop("ROI not yet computed!")
  }
  
  # check input
  to_be_merged <- spe@metadata$roi[[id]]
  
  # remove ROIs first
  if (!is.null(remove.ids)) {
    remove.ids <- as.character(remove.ids)
    to_be_removed <- to_be_merged
    if (sum(!(remove.ids %in% to_be_removed)) > 0.5) {
      stop("Some ROIs are not present in spe@metadata$roi$component. Check remove.ids")
    }
    message(paste("Removing ROIs: ", paste0(remove.ids, collapse = ", "), sep = ""))
    remove_ind <- as.character(spe@metadata$roi[[id]]) %in% remove.ids
    spe@metadata$roi <- spe@metadata$roi[!remove_ind, ]
    to_be_merged <- spe@metadata$roi[[id]]
  }

  if (is.null(merge.list)) {

    message("No ROIs are merged. ")
    ROI_merged <- spe@metadata$roi[[id]]
  
  } else {

    if (any(!(unlist(merge.list) %in% to_be_merged))) {
      stop("Some ROIs are not present in spe@metadata$roi$component or have been removed. 
          Check merge.list and remove.ids!")
    }
    if (any(table(unlist(merge.list)) > 1L)) {
      stop("Each ROI can only be merged once. Check input list!")
    }
    
    list_sizes <- vapply(merge.list, function(rr) {
      length(unique(rr))
    }, numeric(length = 1L))
    if (any(list_sizes < 2L)) {
      stop("Each vector in the list must have at least 2 unique ROI IDs to merge!")
    }
    
    # give names to the list if not input
    if (length(names(merge.list)) < length(merge.list)) {
      names(merge.list) <- vapply(merge.list, function(rr) {
        paste0(sort(rr), collapse = "-")
      }, character(length = 1L))
    }
    
    # merge ROIs
    ROI_merged <- as.character(spe@metadata$roi[[id]])
    for (mm in names(merge.list)) {
      ROI_merged[ROI_merged %in% as.character(merge.list[[mm]])] <- mm
    }
  
  }
  
  if (!rename) {
    # re-order by merged ROI then individual ROIs
    roi_levels <- unique(ROI_merged)
    roi_levels <- c(
      sort(names(merge.list)),
      as.character(sort(
        as.numeric(roi_levels[!(
          roi_levels %in% names(merge.list))])
      ))
    )
    ROI_merged <- factor(ROI_merged, levels = roi_levels)
  } else {
    ROI_merged <- factor(rank(ROI_merged),
                         labels = seq(length(unique(ROI_merged)))
    )
  }
  spe@metadata$roi$component <- ROI_merged
  
  return(spe)
}

