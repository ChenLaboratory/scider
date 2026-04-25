# Perform kernel density estimation on SpatialExperiment for cell types of interest

Perform kernel density estimation on SpatialExperiment for cell types of
interest

## Usage

``` r
gridDensity(
  spe,
  id = if (isVisium != "none") NULL else "cell_type",
  coi = NULL,
  feature = NULL,
  assay = "counts",
  kernel = "gaussian",
  bandwidth = NULL,
  ngrid.x = NULL,
  grid.length.x = NULL,
  diggle = FALSE,
  grid.type = c("hex", "square"),
  isVisium = c("none", "visium", "visiumHD")
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- id:

  A character. The name of the column of colData(spe) containing the
  cell type identifiers. Set to cell_type by default. Set to NULL for
  overall density.

- coi:

  A character vector of cell types of interest (COIs). Default to all
  cell types.

- feature:

  Feature(s) to calculate density with. Must be in rownames(spe).

- assay:

  Name of assay to use for finding feature(s).

- kernel:

  The smoothing kernel. Options are "gaussian", "epanechnikov",
  "quartic" or "disc". For hexagonal grid, only Gaussian is implemented

- bandwidth:

  The smoothing bandwidth. By default performing automatic bandwidth
  selection using cross-validation using function
  spatstat.explore::bw.diggle.

- ngrid.x:

  Number of grids in the x-direction. Ignored when 'grid.length.x' is
  specified. Default to NULL.

- grid.length.x:

  Grid length in the x-direction. If both 'ngrid.x' and 'grid.length.x'
  are NULL, then 'grid.length.x' is set to 100 (micron) by default.

- diggle:

  Logical. If TRUE, use the Jones-Diggle improved edge correction. See
  spatstat.explore::density.ppp() for details.

- grid.type:

  Type of grid can be either hexagon or square.

- isVisium:

  Options of 'none','visium', and 'visiumHD'. If TRUE, converts
  coordinates from pixel to um and fit the density grid to the same
  Visium spots arrangement. visium will use hexagonal while visiumHD
  will use rectangular grid.

## Value

A SpatialExperiment object. Grid density estimates for all cell type of
interest are stored in spe@metadata\$grid_density. Grid information is
stored in spe@metadata\$grid_info

## Examples

``` r
data("xenium_bc_spe")

spe <- gridDensity(spe)
```
