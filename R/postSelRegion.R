#' Merge sel_region from the selectRegion function to SpatialExperiment.
#'
#' @param spe A SpatialExperiment object.
#' @param sel_region A dataframe object. Can be generated from
#' function selectRegion.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#'
#' coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
#'
#' spe <- gridDensity(spe, coi = coi)
#'
#' sel_region <- data.frame(
#'     "node" = seq(10),
#'     "node_x" = seq(10),
#'     "node_y" = seq(10)
#' )
#'
#' spe1 <- postSelRegion(spe, sel_region)
#'
postSelRegion <- function(spe, sel_region) {
    if (!(all(c("node", "node_x", "node_y") %in% colnames(sel_region)))) {
        stop("sel_region should have column node, node_x and node_y.")
    }

    if (!("grid_info" %in% names(spe@metadata))) {
        stop("grid_info is expected in the metadat of spe.")
    }

    components <- S4Vectors::DataFrame(data.frame(
        "component" = "user_defined",
        "members" = sel_region[, "node"],
        "x" = sel_region[, "node_x"],
        "y" = sel_region[, "node_y"]
    ))

    components$xcoord <- spe@metadata$grid_info$xcol[as.numeric(components$x)]
    components$ycoord <- spe@metadata$grid_info$yrow[as.numeric(components$y)]

    spe@metadata$roi <- components

    return(spe)
}
