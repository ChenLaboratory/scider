#' Scale and straighten out Visium coordinates
#' @param spe A SpatialExperiment object.
#' @param distPoint Numeric. Desired point to point distance.
#' @details This function rescale the distance between points to 100um (or other
#' value) to matchthe real distance. In addition, Visium spots have a slight 
#' tilt to them which this function will also fix
#' @export
realignVisium <- function(spe,
                          distPoint = 100) {
  if (is.null(spe$array_col) || is.null(spe$array_row)) {
    stop("spe must have array_col and array_row in its colData.")
  }
  
  # Temporarily scale array_col and array_row by c(1/2,sqrt(3)/2) to match 
  # square arrangement
  old_array_col <- spe@colData$array_col
  old_array_row <- spe@colData$array_row
  spe@colData$array_col <- old_array_col/2
  spe@colData$array_row <- old_array_row/2*sqrt(3)
  spe <- realignVisiumHD(spe,distPoint)
  spe@colData$array_col <- old_array_col
  spe@colData$array_row <- old_array_row
  return(spe)
}

#' Scale and straighten out VisiumHD coordinates
#' @param spe A SpatialExperiment object.
#' @param distPoint Numeric. Desired point to point distance. If NULL, will try 
#' to determine the bin level of spe and use that.
#' @details This function rescale the distance between points to 8um (or other
#' value) to matchthe real distance. In addition, Visium spots have a slight 
#' tilt to them which this function will also fix
#' @export
realignVisiumHD <- function(spe,
                            distPoint = NULL) {
  if (is.null(spe$array_col) || is.null(spe$array_row)) {
    stop("spe must have array_col and array_row in its colData.")
  }
  
  if (is.null(distPoint)) {
    ## Figuring out binsize (2, 8, or 16um).
    distPoint <- .guessVisiumHDBin(spe)
    if (is.null(distPoint)) {
      stop("Couldn't determine distPoint. Please provide one.")
    }
    message(paste0("Using distPoint = ",distPoint))
  }
  
  ## Getting scale factor pxl/um.
  # Use average of 100 distance pairs (one is not too accurate. all is too much).
  n_points <- min(200,ncol(spe)-ncol(spe)%%2)
  coords <- as.matrix(spatialCoords(spe))
  colrow <- as.matrix(spe@colData[,c("array_col","array_row")])
  
  pxl_dist <- coords[1:(n_points/2),]-coords[(n_points/2+1):n_points,]
  pxl_dist <- sqrt(rowSums(pxl_dist^2))
  
  um_dist <- colrow[1:(n_points/2),]-colrow[(n_points/2+1):n_points,]
  um_dist <- sqrt(rowSums(um_dist^2))*distPoint
  scale_factor <- mean(pxl_dist/um_dist)
  
  ## Scaling and realignment
  dim_names <- dimnames(spatialCoords(spe))
  spatialCoords(spe) <- sweep(cbind((colData(spe)$array_col)*distPoint,
                                    (colData(spe)$array_row)*distPoint),
                              2,
                              coords[1,]/scale_factor-colrow[1,]*distPoint,
                              "+")
  dimnames(spatialCoords(spe)) <- dim_names
  
  imgData(spe)$scaleFactor = imgData(spe)$scaleFactor*scale_factor
  return(spe)
}

.guessVisiumHDBin <- function(spe) {
  bin <- NULL
  # check colnames 
  text <- colnames(spe)[1]
  match <- regexpr("^s_(002|008|016)um",text)
  match_loc <- which(match!=-1)[1]
  if (!is.na(match_loc)) {
    bin <- substr(text[match_loc],match[match_loc]+2,match[match_loc]+4)
  } else if (!is.null(spe$bin_size)) {
    # No matches. check for spe$bin_size (VisiumIO).
    bin <- spe$bin_size[1]
  } else {
    # No matches. check colnames of colData (Seurat).
    text <- colnames(spe@colData)
    match <- regexpr("(002|008|016)um",text)
    match_loc <- which(match!=-1)[1]
    if (!is.na(match_loc)) {
      bin <- substr(text[match_loc],match[match_loc],match[match_loc]+2)
    }
  }
  return(as.integer(bin))
}