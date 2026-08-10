# Get top highly variable genes.

Get top highly variable genes.

## Usage

``` r
getHVG(spe, n = 1000, min.total.count = 100, min.prop = 0.01)
```

## Arguments

- spe:

  A SpatialExperiment object.

- n:

  Integer. The number of HVGs.

- min.total.count:

  Numeric. Genes with total counts less than `min.total.count` across
  all cells are not considered for HVGs.

- min.prop:

  Numeric. Genes that have non-zero counts in less than the specified
  proportion (`min.prop`) of all cells are excluded from HVG selection.

## Value

A SpatialExperiment object with HVG information stored in
`rowData(spe)$hvg` as a logical vector. When a mean-dispersion trend is
fitted (i.e. more eligible genes than requested HVGs), the per-gene
trend is also stored in `metadata(spe)$hvg_trend` for use with
`plotHVG`.

## Details

`getHVG` adopts a fast approach of NB dispersion estimation for all the
genes across all cells. A lowess curve is fit to represent
mean-dispersion trend. Top HVGs are selected based on the ratio of
gene-wise dispersion and their trended dispersion.

## Examples

``` r

data("xenium_bc_spe")
spe <- getHVG(spe, n=100)
```
