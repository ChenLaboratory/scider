# Plot reduced dimensions.

plotDR is the main function for plotting reduced dimension. Others are
wrapper functions for convenience.

## Usage

``` r
plotDR(
  spe,
  dimred = NULL,
  dims = c(1, 2),
  group.by = NULL,
  feature = NULL,
  assay = "counts",
  cols = NULL,
  pt.shape = 16,
  pt.size = 1,
  pt.alpha = 0.6,
  label = NULL,
  xlab = NULL,
  ylab = NULL,
  cols.scale = NULL
)

plotUMAP(spe, dimred = "UMAP", ...)

plotPCA(spe, dimred = "PCA", ...)
```

## Arguments

- spe:

  A SpatialExperiment object.

- dimred:

  Name of the reduced dimension in
  [reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)

- dims:

  Numeric vector length 2 for the dimensions to be plotted. Default to
  first two dimensions

- group.by:

  values to group points by. Must be in colData of spe. If NULL, will
  try with 'cols' if available.

- feature:

  Feature to group points by. Must be in rownames(spe).

- assay:

  Name of assay to use for plotting feature.

- cols:

  Colour palette. Can be a vector of colours or a function that accepts
  an integer n and return n colours.

- pt.shape:

  shape of points.

- pt.size:

  size of points.

- pt.alpha:

  alpha of points between 0 and 1.

- label:

  label for the legend

- xlab:

  label for the x-axis

- ylab:

  label for the y-axis

- cols.scale:

  vector of position for color if colors should not be evenly
  positioned. See
  [scale_color_gradientn](https://ggplot2.tidyverse.org/reference/scale_gradient.html).
  Only applicable for continuous values.

- ...:

  Additional arguments pass to plotDR

## Value

A ggplot object.

## Examples

``` r

data("xenium_bc_spe")
spe = runUMAP(spe)
#> PCA not found. Switching to counts assay instead.
plotDR(spe, group.by = "cell_type")

```
