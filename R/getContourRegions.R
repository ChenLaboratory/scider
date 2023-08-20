#' Calculate areas between every two density levels
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#'
#' @return A list of sf objects, each representing the region between two
#' density levels. 
#' @export
#'
#' @examples
#' 
getContourRegions <- function(spe, coi) {
  
  coi_clean <- janitor::make_clean_names(coi)
  coi_clean_contour <- paste(coi_clean, "contour", sep = "_")
  
  contour_data <- spe@metadata[[coi_clean_contour]]
  levs <- unique(contour_data$level)
  nlevs <- length(levs)
  
  # calculate regions for each level
  area_levs <- lapply(1:nlevs, function(ll) contour2sf(spe, 
                                                       contour = coi_clean_contour, 
                                                       coi = coi, 
                                                       level = levs[ll]))
  # area b/w every two levels
  all_areas <- lapply(1:(nlevs - 1), function(ii) 
    st_as_sf(st_union(st_difference(area_levs[[ii]], area_levs[[ii + 1]]))))
  
  # highest level
  all_areas[[nlevs]] <- st_as_sf(area_levs[[nlevs]])
  names(all_areas) <- as.character(1:nlevs)
  
  all_areas
}