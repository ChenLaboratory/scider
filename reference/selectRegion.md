# Select region of interest from plot

Select region of interest from plot

## Usage

``` r
selectRegion(data, x.col = "x", y.col = "y", save.region = FALSE)
```

## Arguments

- data:

  A data.frame object.

- x.col:

  Column name of the x coordinates.

- y.col:

  Column name of the y coordinates.

- save.region:

  Whether to also export the region defined by box/lasso select as an sf
  polygon.

## Value

A data.frame object in the global environment. If save.region is TRUE,
output is a list with a data.frame of selected points and a sf polygon
instead.

## Examples

``` r
data("xenium_bc_spe")

spe_b <- spe[, SummarizedExperiment::colData(spe)$cell_type == "B cells"]

dat <- as.data.frame(SpatialExperiment::spatialCoords(spe_b))

# selectRegion(dat, x.col = "x_centroid", y.col = "y_centroid")
```
