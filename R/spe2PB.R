#' Given a 'SpatialExperiment' data object, create pseudo-bulk
#' samples using the colData information and return a DGEList object
#'
#' @param spe A SpatialExperiment object.
#' @param by.group Logical. Whether to perform pseudo-bulking by group.
#' TRUE by default.
#' @param group.id Character. The column name of the colData(spe) that
#' contains the group information. Default to 'cell_type'.
#' @param by.roi Logical. Whether to perform pseudo-bulking by ROI.
#' TRUE by default.
#' @param roi.only Logical. Whether to filter out pseudo-bulk samples formed
#' by cells not in any ROIs. TRUE by default.
#' @param contour Character. The name of the group or cell type on which
#' the contour level is computed.
#' If NULL, then no pseudo-bulking will be performed based on contour level.
#' Default to NULL.
#'
#' @return An edgeR::DGEList object where each library (column) is a
#' pseudo-bulk sample.
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
#' spe <- findROI(spe, coi = coi)
#'
#' spe <- allocateCells(spe)
#'
#' y <- spe2PB(spe)
#'
spe2PB <- function(spe,
                   by.group = TRUE,
                   group.id = "cell_type",
                   by.roi = TRUE,
                   roi.only = TRUE,
                   contour = NULL) {
    if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
        stop("SpatialExperiment is required but is not installed
         (or can't be loaded)")
    }

    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required but is not installed (or can't be loaded)")
    }

    if (!is(spe, "SpatialExperiment")) {
        stop("spe is not of the SpatialExperiment class")
    }

    # Check 'counts'
    counts <- spe@assays@data$counts
    if (is.null(counts)) stop("spe doesn't contain raw RNA counts")

    # Check 'colData'
    cData <- spe@colData

    grp <- rois <- clvl <- c()

    if (by.group) {
        if (!group.id %in% names(cData)) {
            stop(paste(group.id, "is not found in colData of spe."))
        }
        grp <- cData[, group.id]
    }

    if (by.roi) {
        if (!"roi" %in% names(cData)) {
            message(paste("ROIs are not defined. Proceed without ROIs."))
        } else {
            rois <- paste0("ROI", cData[, "roi"])
        }
    }

    if (!is.null(contour)) {
        contour <- gsub("_contour$", "", contour)
        cont <- paste0(janitor::make_clean_names(contour), "_contour")
        if (!cont %in% names(cData)) {
            stop(paste(
                contour,
                " contour level is not found in colData of spe."
            ))
        } else {
            clvl <- paste0("Lv", cData[, cont])
            if (length(table(clvl)) == 1) {
                warning("Only 1 contour level found. Please check
                whether the contour information is specified correctly.")
            }
        }
    }

    if (is.null(grp) & is.null(rois) & is.null(clvl)) {
        stop("At least one of the following is required:
        group, ROI, or contour information.")
    }

    # Check gene information
    if (ncol(SummarizedExperiment::rowData(spe)) == 0) {
        genes <- data.frame(genes = rownames(spe))
    } else {
        genes <- as.data.frame(SummarizedExperiment::rowData(spe))
    }

    # Pseudo-bulk counts
    group <- paste(grp, rois, clvl, sep = "_")
    group <- factor(gsub("^_|_$", "", group))
    # group_mat <- Matrix::sparse.model.matrix(~ 0 + group)
    group_mat <- stats::model.matrix(~ 0 + group)
    colnames(group_mat) <- gsub("^group", "", colnames(group_mat))
    counts.pb <- counts %*% group_mat

    # Pseudo-bulk sample information
    sample.pb <- data.frame()[seq_len(ncol(counts.pb)), ]
    sample.pb$n.cells <- table(group)

    if (!is.null(grp)) {
        grp.pb <- gsub("_ROI.*$", "", levels(group))
        grp.pb <- gsub("_Lv.*$", "", grp.pb)
        sample.pb$group <- grp.pb
    }
    if (!is.null(rois)) {
        rois.pb <- gsub("^.*ROI", "", levels(group))
        rois.pb <- gsub("_.*$", "", rois.pb)
        sample.pb$ROI <- rois.pb
        if (roi.only) keep <- rois.pb != "None"
    }
    if (!is.null(clvl)) {
        clvl.pb <- gsub("^.*Lv", "", levels(group))
        sample.pb$contour.level <- clvl.pb
    }

    names(sample.pb) <- gsub("group", group.id, names(sample.pb))

    # DGEList
    dge <- edgeR::DGEList(
        counts = as.matrix(counts.pb),
        samples = sample.pb, genes = genes
    )
    if (!is.null(rois) & roi.only) dge <- dge[, keep]

    return(dge)
}
