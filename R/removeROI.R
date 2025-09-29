removeROI <- function (spe, remove.list, id = "component", rename = FALSE) 
{
  if (is.null(spe@metadata$roi)) {
    stop("ROI not yet computed!")
  }
  if (id == "component") {
    to_be_removed <- unique(spe@metadata$roi$component)
  }
  else {
    to_be_removed <- spe@metadata$roi[[id]]
  }
  if (any(!(unlist(remove.list) %in% to_be_removed))) {
    stop("Some ROIs are not present in spe@metadata$roi$component.\n
         Check remove.list")
  }
  if (any(table(unlist(remove.list)) > 1L)) {
    remove.list <- unique(remove.list)
  }
  
  remove_ind <- as.character(spe@metadata$roi$component) %in% as.character(remove.list)
  spe@metadata$roi <- spe@metadata$roi[!remove_ind, ]
  
  if (!rename) {
    
    roi_levels <- levels(spe@metadata$roi$component)
    roi_levels <- setdiff(roi_levels, as.character(remove.list))
    spe@metadata$roi$component <- factor(spe@metadata$roi$component, levels = roi_levels)
    
  }
  else {
    
    spe@metadata$roi$component <- factor(rank(spe@metadata$roi$component), 
                                         labels = seq(length(unique(spe@metadata$roi$component))))
    
  }
  
  return(spe)
  
}
