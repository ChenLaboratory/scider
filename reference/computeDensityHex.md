# Perform kernel density estimation on SpatialExperiment

Perform kernel density estimation on SpatialExperiment

## Usage

``` r
computeDensityHex(
  xy,
  kernel = c("gaussian"),
  bandwidth = NULL,
  weights = NULL,
  ngrid.x = NULL,
  xlim = NULL,
  ylim = NULL,
  diggle = FALSE,
  gridInfo = FALSE
)
```

## Arguments

- xy:

  A numeric matrix of spatial coordinates.

- kernel:

  The smoothing kernel. Options are gaussian, epanechnikov, quartic or
  disc. ONLY GAUSSIAN IS IMPLEMENTED

- bandwidth:

  The smoothing bandwidth. By default performing automatic bandwidth
  selection using cross-validation using function
  spatstat.explore::bw.diggle.

- weights:

  Optional weights to be attached to the points.

- ngrid.x:

  Number of grids in the x-direction.

- xlim:

  The range of the x-coordinates of the image.

- ylim:

  The range of the y-coordinates of the image.

- diggle:

  Logical. If TRUE, use the Jones-Diggle improved edge correction. See
  spatstat.explore::density.ppp() for details.

- gridInfo:

  Logical. If TRUE, then the grid information is also returned.

## Value

Output from spatstat.explore::density.ppp.
