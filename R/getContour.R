#' Get contour from density
#'
#' @param spe A SpatialExperiment object.
#' @param coi A character vector of cell types of interest (COIs).
#' All cell types are chosen if NULL or `overall`.
#' @param equal.cell Logical. Whether to use produce contour levels so that
#' there are roughly the same number of cells of the COI at each level. 
#' Default to TRUE.
#' @param bins An integer. Number of contour levels.
#' @param binwidth A numeric scale of the smoothing bandwidth.
#' @param breaks A numeric scale referring to the breaks in
#' `ggplot2:::contour_breaks`.
#' @param id A character. The name of the column of colData(spe) containing
#' the cell type identifiers. Set to cell_type by default or in_tissue if spe 
#' is Visium. Only needed when \code{equal.cell = TRUE}. 
#'
#' @return A SpatialExperiment object. An sf object of the contour region of
#' the specified level is stored in the metadata of the
#' SpatialExperiment object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' spe <- gridDensity(spe)
#'
#' coi <- "Breast cancer"
#'
#' spe <- getContour(spe, coi = coi)
#'
getContour <- function(spe, coi = NULL, equal.cell = TRUE, bins = NULL,
                       binwidth = NULL, breaks = NULL, 
                       id = NULL) {
    
    if (is.null(spe@metadata$grid_density)) {
        stop("Have to calculate grid density, run gridDensity() first!")
    }
  
    if (is.null(id)) {
      id = `if`(!is.null(spe@metadata$grid_info$isVisium),
                "in_tissue",
                "cell_type")
    }
  
    if (equal.cell && !id %in% colnames(colData(spe))) {
        stop(paste(id, "is not a column of the colData."))
    }

    if ( !is.null(coi) & !("overall" %in% coi) ){
        if (length(which(!coi %in% names(table(colData(spe)[[id]])))) > 0L) {
            stop(paste(paste0(
                coi[which(!coi %in%
                    names(table(colData(spe)[[id]])))],
                collapse = ", "
            ), "not found in data!", sep = " "))
        }
    } else coi <- "overall"

    coi_clean <- janitor::make_clean_names(coi)
    dens_cols <- paste("density", coi_clean, sep = "_")

    # grid level density data
    dens <- spe@metadata$grid_density
    dups <- duplicated(dens[, c("y_grid", "x_grid"), drop = FALSE],
        fromLast = TRUE
    )
    dens <- dens[!dups, , drop = FALSE]
    dens <- as.data.frame(dens)

    if (!all(dens_cols %in% colnames(dens))) {
        stop("Density of COI is not yet computed.")
    }

    if (length(dens_cols) > 1L) {
        message("Finding contour using total density of input COIs. ")
        dens$density_coi <- rowSums(dens[, which(colnames(dens) %in%
            dens_cols),
        drop = FALSE
        ])
    } else {
        dens$density_coi <- dens[, dens_cols]
    }
    
    dens <- dens[, c(seq_len(5), which(colnames(dens) ==
        "density_coi"))]

    # filter out negative densities when calculating contours
    dens <- dens[dens$density_coi > 0L, ]

    # levels for contour
    if (!equal.cell) {
        if (is.null(bins) && is.null(binwidth) && is.null(breaks)) {
        message("Using bins = 10 to draw contours.")
        bins <- 10L
        }
        if (!is.null(bins)) binwidth <- breaks <- NULL
        if (is.null(bins) && !is.null(binwidth)) breaks <- NULL
    } else {
        # calculate density level breaks to get roughly the same number of cells at each level
        if (is.null(bins)) {
            message("Using bins = 10 to draw contours with equal cell numbers.")
            bins <- 10L
            binwidth <- breaks <- NULL
        } else {
            message(paste("Using bins =", bins, "to draw contours with equal cell numbers.", sep = " "))
        }
        ## count no of cells of coi in each grid
        ## note this no can be very different from the expected no
        coi_coords <- as.data.frame(spatialCoords(spe))
        if(!"overall" %in% coi){
            coi_coords <- coi_coords[colData(spe)[[id]] %in% coi, ]
        }

        if (spe@metadata$grid_info$grid_type == "hex") {
          hcellsID = hexDensity::xy2hcell(x=coi_coords$x_centroid,y=coi_coords$y_centroid,
                                          xbins=spe@metadata$grid_info$xbins,
                                          xbnds=spe@metadata$grid_info$xlim,
                                          ybnds=spe@metadata$grid_info$ylim,
                                          shape=spe@metadata$grid_info$shape)
          coi_coords$hcellsID=hcellsID
          coi_coords$x_node = (hcellsID-1)%%spe@metadata$grid_info$dims[1]+1
          coi_coords$y_node = (hcellsID-1)%/%spe@metadata$grid_info$dims[1]+1
        } else {
          coi_coords$x_node <- vapply(coi_coords$x_centroid, function(xx) {
              which.min(abs(spe@metadata$grid_info$xcol - xx))
          }, numeric(1))
          coi_coords$y_node <- vapply(coi_coords$y_centroid, function(yy) {
              which.min(abs(spe@metadata$grid_info$yrow - yy))
          }, numeric(1))
        }
        coi_coords$node <- paste(coi_coords$x_node, coi_coords$y_node, sep = "-")
        freq <- c(table(coi_coords$node))
        dens$n <- freq[dens$node]
        dens$n <- ifelse(is.na(dens$n), 0L, dens$n)
        dens_expanded <- rep(dens$density_coi, times = dens$n)
        dens_expanded <- dens_expanded[dens_expanded > 0L]
        qq <- seq(0, 1, length.out = bins + 1)[-1]
        if (qq[length(qq)] == 1L) qq <- qq[-length(qq)]
        breaks <- as.vector(quantile(dens_expanded, probs = qq))
        binwidth <- bins <- NULL
    }
    

    # note that when calculating contours, density is not filtered at any
    # quantile cutoff!

    contour <- compute_group0(dens,
        z.range = range(dens$density_coi, na.rm = TRUE, finite = TRUE),
        bins = bins,
        binwidth = binwidth,
        breaks = breaks,
        na.rm = FALSE,
        grid_type = spe@metadata$grid_info$grid_type
    )

    contour$level <- as.factor(as.numeric(as.factor(contour$cutoff)))

    coi_clean_output <- ifelse(length(coi_clean) == 1L, coi_clean, paste(sort(coi_clean), collapse="_"))
    spe@metadata[[paste(coi_clean_output,
        "contour",
        sep = "_"
    )]] <- S4Vectors::DataFrame(contour)

    return(spe)
}



#### lower level functions for computing the contours.
# compute contour groups (grabbed from ggplot2)
xyz_to_isolines_square <- function(data, breaks) {
  # Convert vector of data to raster z
  x_pos <- as.integer(factor(data$x_grid,
                             levels =
                               sort(unique00(data$x_grid))
  ))
  y_pos <- as.integer(factor(data$y_grid,
                             levels =
                               sort(unique00(data$y_grid))
  ))
  nrow <- max(y_pos)
  ncol <- max(x_pos)
  z <- matrix(NA_real_, nrow = nrow, ncol = ncol)
  z[cbind(y_pos, x_pos)] <- data$density_coi
  
    isoband::isolines(
        x = sort(unique00(data$x_grid)),
        y = sort(unique00(data$y_grid)),
        z = z,
        levels = breaks
    )
}
xyz_to_isolines_hex = function(data, breaks) {
  x.coords=sort(unique00(data$x_grid))
  x.coords.left=x.coords[seq.int(1,length(x.coords),2)]
  x.coords.right=x.coords[seq.int(2,length(x.coords),2)]
  y.coords=sort(unique00(data$y_grid))
  # Convert vector of data to raster
  ncol=diff(range(data$node_x))+1
  nrow=diff(range(data$node_y))+1
  z = matrix(NA_real_, nrow = nrow, ncol = ncol)
  z[cbind(data$node_y, data$node_x)] <- data$density_coi
  isolines=hexDensity::meanderingTriangles(x.coords.left,x.coords.right,
                                           y.coords,z,breaks)

  return(isolines)
}

iso_to_path <- function(iso, group = 1) {
    lengths <- vapply(iso, function(x) length(x$x), integer(1))

    if (all(lengths == 0)) {
        message("Zero contours were generated.")
        return(data_frame00())
    }

    levels <- names(iso)
    xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
    ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
    ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
    item_id <- rep(seq_along(iso), lengths)

    # Add leading zeros so that groups can be properly sorted
    groups <- paste(group, sprintf("%03d", item_id), sprintf("%03d", ids),
        sep = "-"
    )
    groups <- factor(groups)

    data_frame00(
        level = rep(levels, lengths),
        x = xs,
        y = ys,
        piece = as.integer(groups),
        group = groups,
        .size = length(xs)
    )
}

compute_group0 <- function(data, z.range, bins = NULL,
                          binwidth = NULL, breaks = NULL, na.rm = FALSE,
                          grid_type) {
    breaks <- contour_brks(z.range, bins, binwidth, breaks)
    isolines <- `if` (grid_type=="hex",
                      xyz_to_isolines_hex(data, breaks),
                      xyz_to_isolines_square(data, breaks))
    path_df <- iso_to_path(isolines, data$group[1])
    path_df$cutoff <- as.numeric(path_df$level)
    # path_df$level <- as.numeric(path_df$level)
    # path_df$nlevel <- scales::rescale_max(path_df$level)
    path_df
}

