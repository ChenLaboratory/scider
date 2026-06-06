# Perform log normalization for counts

Perform log normalization for counts

## Usage

``` r
normalizeAssay(
  spe,
  transformation = c("log"),
  scale.factor = 10000,
  assay = "counts",
  name = "logcounts"
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- transformation:

  Choice of transformation. "Log" for log1p

- scale.factor:

  Factor to multiply the count of each cell by. A single value or a
  numeric vector of length equal to the number of cells.

- assay:

  Name of assay in spe to perform the transformation on

- name:

  Name of the transformed assay

## Value

A SpatialExperiment object

## Examples

``` r
data("xenium_bc_spe")
spe <- normalizeAssay(spe)
#> Error in `assays<-`(`*tmp*`, withDimnames = withDimnames, ..., value = `*vtmp*`): please use 'assay(x, withDimnames=FALSE) <- value' or 'assays(x,
#>   withDimnames=FALSE) <- value' when the rownames or colnames of the supplied
#>   assay(s) are not identical to those of the receiving SpatialExperiment object
#>   'x'
```
