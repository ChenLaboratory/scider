# Construct a distance-based neighbour list from cell coordinates.

Construct a distance-based neighbour list from cell coordinates.

## Usage

``` r
findNbrsSpatial(
  spe,
  k = NULL,
  radius = NULL,
  dist_func = c("idw", "exp", "binary", "none"),
  standardisation = c("none", "row"),
  scale = 1,
  nbrs_name = NULL,
  cpu_threads = 6
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- k:

  Integer scalar for number of nearest neighbours to find. Can be used
  with radius. See details.

- radius:

  Numeric for maximum distance to search for neighbours. Can be with k.
  See details

- dist_func:

  Options for distance-based weight. "idw" for inverse distance, "exp"
  for exponential decay, "binary" for constant weight, and "none" for
  raw euclidean distance.

- standardisation:

  Options for weight standardisation. "none" for nothing, and "row" for
  dividing weights by number of neighbours.

- scale:

  Numeric scaler for weight scaling.

- nbrs_name:

  Name of the neighbour list to be stored. Default to be "spatial".

- cpu_threads:

  Number of cpu threads for parallel computation.

## Value

A SpatialExperiment object with neighbour list stored in
`spe@metadata$nbrs$cell[[nbrs_name]]`

## Details

if only `k` is provided, neighbours are found using
[findKNN](https://rdrr.io/pkg/BiocNeighbors/man/findKNN.html). If only
`radius` is provided, neighbours are found using
[findNeighbors](https://rdrr.io/pkg/BiocNeighbors/man/findNeighbors.html).
If both are provided, then knn is done first then neighbours are
filtered to only those within radius.

## Examples

``` r

data("xenium_bc_spe")
spe <- findNbrsSpatial(spe,k=20,radius=100)
```
