# Construct a neighbour list from grid coordinates.

Construct a neighbour list from grid coordinates.

## Usage

``` r
findNbrsGrid(
  spe,
  n = 1,
  radius = NULL,
  diagonal = FALSE,
  dist_func = c("idw", "exp", "binary", "none"),
  dist_type = c("euclidean", "manhattan"),
  standardisation = c("row", "none"),
  scale = 1,
  nbrs_name = NULL,
  cpu_threads = 6
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- n:

  Integer. Search for neighbours within (...). Either the number of
  neighbors or radius

- radius:

  Numeric. Search for neighbours within the radius.

- diagonal:

  Whether to consider diagonal connection if using square grid

- dist_func:

  Options for distance-based weight. "idw" for inverse distance, "exp"
  for exponential decay, "binary" for constant weight, and "raw" for raw
  distance.

- dist_type:

  Options of using euclidean or manhattan for distance calculation

- standardisation:

  Options for weight standardisation. "none" for nothing, and "row" for
  dividing weights by number of neighbours.

- scale:

  Numeric scaler for weight scaling.

- nbrs_name:

  Name of the neighbour list to be stored. Default to be "grid".

- cpu_threads:

  Number of cpu threads for parallel computation.

## Value

A SpatialExperiment object with neighbour list stored in
`spe@metadata$nbrs$grid[[nbrs_name]]`

## Details

If n is used, distance is scaled to unit distance

## Examples

``` r
data("xenium_bc_spe")
spe <- gridDensity(spe)
spe <- findNbrsGrid(spe,n=3)
```
