# Find ROIs based on cell type-specific densities via graph-based method.

Find ROIs based on cell type-specific densities via graph-based method.

## Usage

``` r
findROI(
  spe,
  coi = NULL,
  probs = 0.85,
  min.density = NULL,
  ngrid.min = 20,
  method = c("greedy", "walktrap", "connected", "hdbscan", "eigen", "dbscan"),
  diag.nodes = FALSE,
  sequential.roi.name = TRUE,
  zoom.in = FALSE,
  zoom.in.size = 500L,
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- coi:

  A character vector of cell types of interest (COIs). Default to all
  cell types.

- probs:

  A numeric scalar. The threshold of proportion that used to filter
  grids by density. Default to 0.85.

- min.density:

  A numeric value. The cut-off value used to filter grids by density.
  Default is NULL and overwrites probs.

- ngrid.min:

  An integer. The minimum number of grids required for defining a ROI.
  Default to 20.

- method:

  The community dectection method to be used, possible options are
  greedy, walktrap, connected, hdbscan, eigen or dbscan. Default to
  greedy, can be abbreviated.

- diag.nodes:

  Logical. Set this to TRUE to allow diagonal grid points to be adjacent
  nodes.

- sequential.roi.name:

  Logical. Set this to FALSE if you want the original ROI name before
  filtering are retained.

- zoom.in:

  Logical. For very large ROIs, whether to zoom in and try to get more
  refined ROIs.

- zoom.in.size:

  A numeric scaler. Smallest size of an ROI to be able to zoom in.
  Default is 500L.

- ...:

  Other parameters that passed to walktrap.community when method =
  "walktrap".

## Value

A SpatialExperiment object.

## Examples

``` r
data("xenium_bc_spe")

coi <- c("Breast cancer", "Fibroblasts")

spe <- gridDensity(spe, coi = coi)

spe <- findROI(spe, coi = coi, method = "walktrap")
```
