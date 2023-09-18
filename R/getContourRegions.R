#' Calculate areas between every two density levels
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#'
#' @return A list of sf objects, each representing the region between
#' two contour density levels.
#'
getContourRegions <- function(spe, coi) {
  coi_clean <- janitor::make_clean_names(coi)
  coi_clean_contour <- paste(coi_clean, "contour", sep = "_")

  contour_data <- spe@metadata[[coi_clean_contour]]
  levs <- unique(contour_data$cutoff)
  nlevs <- length(levs)

  # calculate regions for each level
  area_levs <- lapply(seq_len(nlevs), function(ll) {
    contour2sf(spe,
      contour = coi_clean_contour,
      coi = coi,
      cutoff = levs[ll]
    )
  })
  # area b/w every two levels
  all_areas <- lapply(seq_len((nlevs - 1)), function(ii) {
    sf::st_as_sf(sf::st_union(sf::st_difference(
      area_levs[[ii]],
      area_levs[[ii + 1]]
    )))
  })

  # highest level
  all_areas[[nlevs]] <- sf::st_as_sf(area_levs[[nlevs]])
  names(all_areas) <- as.character(seq_len(nlevs))

  all_areas
}
