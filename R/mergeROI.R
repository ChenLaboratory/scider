#' Manually merge ROIs
#'
#' @param spe A SpatialExperiment object.
#' @param merge.list A (named) list of vectors of ROI ids to be merged. Each
#' vector in the list should be of length greater than or equal to 2. If no
#' name is specified, the merged ROI will be named by concatenating ROIs being
#' merged. 
#' @param id Character. The name of the column in \code{spe@metadata$roi} that
#' stores the ROIs to be merged. Default is "component". 
#' @param rename Logical. If TRUE, names of merge.list are ignored. ROIs will be
#' given a new name. For the unmerged ROIs, their new names are not necessarily
#' the same as those before merging. 
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
mergeROI <- function(spe, merge.list, id = "component", rename = FALSE) {
  
  # check ROI exists
  if (is.null(spe@metadata$roi)) {
    stop("ROI not yet computed!")
  }
  
  # check input
  if (id == "component") {
    to_be_merged <- unique(spe@metadata$roi$component)
  } else {
    to_be_merged <- spe@metadata$roi[[id]]
  }
  if (any(!(unlist(merge.list) %in% to_be_merged))) {
    stop("Some ROIs are not present in spe@metadata$roi$component. Check input list!")
  }
  if (any(table(unlist(merge.list)) > 1L)) {
    stop("Each ROI can only be merged once. Check input list!")
  }
  list_sizes <- vapply(merge.list, function (rr) {
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
  
  # keep a copy of the original ROI list
  if (id == "component") {
    ROI_merged <- spe@metadata$roi[["component_before_merge"]] <- spe@metadata$roi$component
  } else {
    ROI_merged <- spe@metadata$roi[[id]]
  }
  
  # merge ROIs
  ROI_merged <- as.character(ROI_merged)
  for (mm in names(merge.list)) {
    ROI_merged[ROI_merged %in% as.character(merge.list[[mm]])] <- mm
  }
  if (!rename) {
    # re-order by merged ROI then individual ROIs
    roi_levels <- unique(ROI_merged)
    roi_levels <- c(sort(names(merge.list)), 
                    as.character(sort(as.numeric(roi_levels[!(roi_levels %in% names(merge.list))]))))
    ROI_merged <- factor(ROI_merged, levels = roi_levels)
    spe@metadata$roi$component <- ROI_merged
  } else {
    ROI_merged <- factor(rank(ROI_merged),
                         labels = seq(length(unique(ROI_merged))))
    spe@metadata$roi$component <- ROI_merged
  }
  
  return(spe)
  
}




