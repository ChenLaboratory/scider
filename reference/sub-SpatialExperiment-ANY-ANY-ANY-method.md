# Subset for grid level analysis

Overwrite the default SpatialExperiment subsetting method to ensure
'grid_density' is also subsetted if 'gridLevelAnalysis' is TRUE (1
polygon 1 spot)

## Usage

``` r
# S4 method for class 'SpatialExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A SpatialExperiment object.

- i:

  row indices for subsetting.

- j:

  col indices for subsetting.

- ...:

  further arguments to be passed to or from other methods.

- drop:

  passed on to \[ indexing operator.

## Value

A SpatialExperiment object.
