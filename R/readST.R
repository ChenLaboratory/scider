# Parse a 10x scalefactors_json.json file without a JSON dependency.
# The file is a flat, single-level object of numeric key:value pairs, so a
# small regex-based parser is sufficient. Returns a named list of numerics.
.readScaleFactors <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "")
  txt <- gsub("[{}\"[:space:]]", "", txt)
  pairs <- strsplit(strsplit(txt, ",")[[1]], ":")
  vals <- as.numeric(vapply(pairs, `[`, "", 2))
  names(vals) <- vapply(pairs, `[`, "", 1)
  as.list(vals)
}

#' Read Visium output into spe
#' @param dir directory containing the Visium files
#' @param sample_id Name of the sample.
#' @param count Name of the h5 file with the count assay.
#' @param coord Path to the tissue coordinates file (csv or parquet), or a
#' data.frame of coordinates (rownames = barcodes) with columns
#' 'pxl_col_in_fullres' and 'pxl_row_in_fullres'.
#' @param image Names of the image files.
#' @param scale_factors Names of the scale factors file
#' @param feature_type Feature type to retain. Defaults to "Gene Expression" to
#' exclude non-gene features. Set to NULL to keep all features.
#' @param pixel_to_micron Logical. If TRUE (default), convert the spot
#' coordinates from full-resolution pixels to microns using
#' 'spot_diameter_fullres' from the scale factors file, and store the
#' conversion factor in metadata(spe)$um_per_pixel. Set to FALSE to keep the
#' coordinates in pixels (previous behaviour).
#' @export
readVisium <- function(dir,
                       sample_id="sample01",
                       count = NULL,
                       coord = NULL,
                       image = NULL,
                       scale_factors = NULL,
                       feature_type = "Gene Expression",
                       pixel_to_micron = TRUE) {
  # read in counts
  count <- count %||% file.path(dir,"filtered_feature_bc_matrix.h5")
  sce <- DropletUtils::read10xCounts(count, col.names = TRUE)

  # Filter to desired feature type and set Symbol as rownames
  rd <- SummarizedExperiment::rowData(sce)
  if (!is.null(feature_type)) {
    keep <- rd$Type %in% feature_type
    sce <- sce[keep, ]
    rd <- rd[keep, ]
  }

  # Materialise into an in-memory sparse matrix so the SPE is self-contained
  # and can be saved/loaded without the original .h5 file being present.
  counts_mat <- methods::as(SummarizedExperiment::assay(sce), "dgCMatrix")
  rownames(counts_mat) <- rd$Symbol
  rownames(rd) <- rd$Symbol

  # Per-cell QC metrics computed directly from the count matrix.
  n_counts <- Matrix::colSums(counts_mat)
  n_genes  <- stats::setNames(as.integer(diff(counts_mat@p)), colnames(counts_mat))

  # read in image
  image <- image %||% file.path(dir, "spatial",
                                c("tissue_lowres_image.png",
                                  "tissue_hires_image.png"))
  scale_factors <- scale_factors %||% file.path(dir,"spatial","scalefactors_json.json")
  sf_json <- .readScaleFactors(scale_factors)
  # Prefer an explicit microns-per-pixel (present in VisiumHD scale factors);
  # otherwise derive it from the 55 um spot diameter (standard Visium).
  um_per_pixel <- if (!pixel_to_micron) {
    NULL
  } else if (!is.null(sf_json$microns_per_pixel)) {
    sf_json$microns_per_pixel
  } else if (!is.null(sf_json$spot_diameter_fullres)) {
    55 / sf_json$spot_diameter_fullres
  } else NULL
  # load=TRUE embeds the image pixels in the SPE object so the RDS is
  # self-contained and portable (load=FALSE only stores the file path, which
  # breaks when the RDS is moved to a machine without access to the original
  # directory).
  img <- SpatialExperiment::readImgData(
    path=dir,
    imageSources = image,
    scaleFactors = scale_factors,
    sample_id=sample_id,
    load=TRUE)

  # read in coords. 'coord' may be a file path (csv/parquet) or a pre-built
  # data.frame of coordinates (rownames = barcodes, with pxl_col_in_fullres /
  # pxl_row_in_fullres columns), as used by the segmented VisiumHD path.
  coord <- coord %||% file.path(dir,"spatial","tissue_positions.csv")
  if (is.data.frame(coord)) {
    spatial <- coord
  } else if (grepl(".csv$",coord)) {
    spatial <- read.csv(coord,row.names=1)
  } else if (grepl(".parquet$",coord)) {
    spatial <- as.data.frame(arrow::read_parquet(coord))
    rownames(spatial) = spatial$barcode
  }
  matches <- intersect(colnames(sce), rownames(spatial))
  spatial <- spatial[matches, ]
  # Keep the count matrix in step with the matched barcodes/cells.
  counts_mat <- counts_mat[, matches, drop = FALSE]
  # Attach QC metrics, matching on barcode to guard against row reordering.
  spatial$n_counts <- as.integer(n_counts[rownames(spatial)])
  spatial$n_genes  <- as.integer(n_genes[rownames(spatial)])

  if (!is.null(um_per_pixel)) {
    spatial$pxl_col_in_fullres <- spatial$pxl_col_in_fullres * um_per_pixel
    spatial$pxl_row_in_fullres <- spatial$pxl_row_in_fullres * um_per_pixel
  }

  spe <- SpatialExperiment(
    assays = list(counts = counts_mat),
    rowData = rd,
    colData = S4Vectors::DataFrame(spatial),
    spatialCoordsNames = c("pxl_col_in_fullres","pxl_row_in_fullres"),
    imgData=img,
    sample_id=sample_id
  )
  # Always record the coordinate unit so downstream code never has to guess.
  # um_per_pixel is only stored when a conversion was actually applied.
  if (!is.null(um_per_pixel)) {
    spe@metadata$um_per_pixel <- um_per_pixel
    spe@metadata$coord_unit <- "micron"
  } else {
    spe@metadata$coord_unit <- "pixel"
  }
  return(spe)
}

# Build a per-cell coordinate table from a VisiumHD cell-segmentation GeoJSON,
# using each cell polygon's centroid as its full-resolution pixel coordinate.
# Cell IDs are mapped to count-matrix barcodes as 'cellid_<zero-padded 9>-1',
# matching Space Ranger's segmented output.
.readSegmentedCoords <- function(geojson) {
  gdf  <- sf::st_read(geojson, quiet = TRUE)
  # The polygons are in planar full-resolution pixel coordinates, but GeoJSON is
  # WGS84 by spec, so sf tags them EPSG:4326 and routes centroids through the
  # spherical s2 engine, which mis-reads pixels as lon/lat and rejects valid
  # cells ("Loop is not valid"). Drop the CRS so planar (GEOS) geometry is used.
  geom <- sf::st_set_crs(sf::st_geometry(gdf), NA)
  cent <- sf::st_coordinates(suppressWarnings(sf::st_centroid(geom)))
  barcode <- sprintf("cellid_%09d-1", as.integer(gdf$cell_id))
  data.frame(
    pxl_col_in_fullres = cent[, 1],
    pxl_row_in_fullres = cent[, 2],
    row.names          = barcode
  )
}

#' Read VisiumHD output into spe
#' @param dir directory containing the VisiumHD files
#' @param bin Which output to read. Bin sizes "016um", "008um", "002um" read the
#' corresponding 'binned_outputs/square_*' folder. "segmented" reads the
#' cell-segmentation results in 'segmented_outputs': the count matrix
#' 'filtered_feature_cell_matrix.h5', with per-cell coordinates taken from the
#' centroids of 'cell_segmentations.geojson'.
#' @param ... Parameters for readVisium
#' @details For "segmented", cells have no array_row/array_col grid, so the
#' result is cell-level (like Xenium) and is not compatible with the Visium
#' spot-grid options of gridDensity() / trimEdge().
#' @export
readVisiumHD <- function(dir,
                         bin = c("016um", "008um", "002um", "segmented"),
                         ...) {
  bin <- match.arg(bin)
  args <- list(...)

  if (bin == "segmented") {
    seg_dir    <- file.path(dir, "segmented_outputs")
    args$dir   <- seg_dir
    args$count <- args$count %||%
      file.path(seg_dir, "filtered_feature_cell_matrix.h5")
    if (is.null(args$coord)) {
      args$coord <- .readSegmentedCoords(
        file.path(seg_dir, "cell_segmentations.geojson"))
    }
  } else {
    bin_dir    <- file.path(dir, "binned_outputs", paste0("square_", bin))
    args$dir   <- bin_dir
    args$coord <- args$coord %||%
      file.path(bin_dir, "spatial", "tissue_positions.parquet")
  }

  return(do.call(readVisium, args))
}

#' Read Xenium output into spe
#' @param dir directory containing the Xenium files
#' @param sample_id Name of the sample.
#' @param count Name of the h5 file with the count assay.
#' @param coord Name of the parquet file with the tissue coordinates
#' @param image Names of the ome.tif image files.
#' @param image_reso resolution of the image to use. From 1-8 (lower = better resolution).
#' See https://kb.10xgenomics.com/hc/articles/11636252598925. Default to 6
#' @param image_layer Which layer of the tiff image to use. Default is the
#' middle-most layer
#' @param feature_type Feature type to retain. Defaults to "Gene Expression" to
#' exclude control codewords. Set to NULL to keep all features.
#' @export
readXenium <- function(dir,
                       sample_id="sample01",
                       count = file.path(dir,"cell_feature_matrix.h5"),
                       coord = file.path(dir,"cells.parquet"),
                       image = file.path(dir,"morphology.ome.tif"),
                       image_reso = 6,
                       image_layer = NULL,
                       feature_type = "Gene Expression") {
  ## read in count
  sce <- DropletUtils::read10xCounts(count, col.names = TRUE)

  ## Filter to desired feature type and set Symbol as rownames
  rd <- SummarizedExperiment::rowData(sce)
  if (!is.null(feature_type)) {
    keep <- rd$Type %in% feature_type
    sce <- sce[keep, ]
    rd <- rd[keep, ]
  }

  # Materialise into an in-memory sparse matrix so the SPE is self-contained
  # and can be saved/loaded without the original .h5 file being present.
  counts_mat <- methods::as(SummarizedExperiment::assay(sce), "dgCMatrix")
  rownames(counts_mat) <- rd$Symbol
  rownames(rd)         <- rd$Symbol

  # Per-cell QC metrics computed directly from the count matrix.
  n_counts <- Matrix::colSums(counts_mat)
  n_genes  <- stats::setNames(as.integer(diff(counts_mat@p)), colnames(counts_mat))

  ## read in coords.
  spatial <- as.data.frame(arrow::read_parquet(coord))
  rownames(spatial) <- spatial$cell_id
  spatial <- spatial[colnames(counts_mat), ]
  spatial$n_counts <- as.integer(n_counts)
  spatial$n_genes  <- as.integer(n_genes)
  
  ## read in image
  imgData <- tryCatch({
    image_reso <- min(image_reso,length(RBioFormats::read.metadata(image)))
    img <- RBioFormats::read.image(image,resolution=image_reso,
                                   read.metadata = FALSE,
                                   normalize = FALSE,
                                   proprietary.metadata = FALSE)
    # Use the middle layer by default
    i <- image_layer %||% dim(img@.Data)[3]%/%2
    # Normalize pixel value to max of 1
    img_raster <- t(img@.Data[,,i])
    img_raster <- grDevices::as.raster(img_raster/max(img_raster))
    # Map ome.tiff resolution to scaling factor. 
    # See https://kb.10xgenomics.com/hc/articles/11636252598925
    scaling <- 1/c(0.2125,0.4250,0.8500,1.7000,3.4000,6.8000,13.6000,27.2000)
    S4Vectors::DataFrame(sample_id=sample_id,
                                image_id="temp",
                                data=I(list(SpatialImage(img_raster))),
                                scaleFactor=scaling[image_reso])
  }, error = function(e) {
    if (grepl("java.lang",e$message,fixed=TRUE)) {
      cat("Skipping image due to Java out of memory error. Try increasing", 
          "maximum heap size before importing scider. (e.g.: ",
          "options(java.parameters = '-Xmx4g') to increase heap size to 4gb)")
    }
    return(NULL)
  }
  )
  
  
  spe <- SpatialExperiment(
    assays = list(counts = counts_mat),
    rowData = rd,
    colData = S4Vectors::DataFrame(spatial),
    spatialCoordsNames = c("x_centroid","y_centroid"),
    imgData=imgData,
    sample_id=sample_id
  )
}

########## Proseg ##########
#' Read Proseg V2 output into spe
#' @param dir directory containing the Proseg files
#' @param sample_id Name of the sample.
#' @param count Name of the file with the count assay.
#' @param coord Name of the file with the tissue coordinates
#' @param coordNames Name of the coordinates for the spe
#' @param gene Name of the file with the gene metadata. This 
#' @export
#' @details This does not work on zarr file output of proseg V3
readProseg <- function(dir,
                       sample_id="sample01",
                       count = "expected-counts\\.(csv|mtx)\\.gz",
                       coord = "cell-metadata.csv.gz",
                       gene = "gene-metadata.csv.gz",
                       coordNames = c("centroid_x","centroid_y","centroid_z")) {
  # Check if required files exist
  coord <- file.path(dir,coord)
  if (!file.exists(coord)) stop(paste0("Couldn't find ",coord))
  count <- list.files(path=dir,pattern=count,full.names = TRUE)[1]
  if (length(count)==0) stop("Couldn't find a suitable expected-counts file.")
  
  # Read in gene metadata. Do this before count 
  gene_file <- file.path(dir,gene)
  if (file.exists(gene_file)) {
    gene <- S4Vectors::DataFrame(read.csv(gzfile(gene_file)))
  } else {
    if (!missing(gene)) message("Couldn't find the file containing gene metadata. Skipping")
    gene <- NULL
  }
  
  # read in counts
  switch(sub(".*\\.(csv\\.gz|mtx\\.gz)$", "\\1",count),
         "csv.gz" = {
           count <- read.csv(gzfile(count))
           },
         "mtx.gz" = {
           count <- Matrix::readMM(gzfile(count))
           if (is.null(gene)) warning("mtx file does not have dimnames.")
           else colnames(count) <- gene$gene
           },
         stop("count need to be either a 'csv.gz' or a 'mtx.gz' file"))

  # Build genes x cells sparse matrix and compute per-cell QC metrics.
  counts_mat <- methods::as(Matrix::t(count), "dgCMatrix")
  n_counts <- as.integer(Matrix::colSums(counts_mat))
  n_genes  <- as.integer(diff(counts_mat@p))

  # read in coords
  coord <- S4Vectors::DataFrame(read.csv(gzfile(coord)))
  if (!(is.null(coord$cluster) || "cluster"%in% coordNames)) {
    coord$cluster <- as.factor(coord$cluster)
  }
  coord$n_counts <- n_counts
  coord$n_genes  <- n_genes

  spe <- SpatialExperiment(
    assays = list(counts = counts_mat),
    colData = coord,
    rowData = gene,
    spatialCoordsNames = coordNames,
    sample_id=sample_id
  )
  return(spe)
}