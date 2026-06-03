# Plot cells based on spatial coordinates.

Plot cells based on spatial coordinates.

## Usage

``` r
plotSpatial(
  spe,
  group.by = NULL,
  feature = NULL,
  assay = "counts",
  type = c("raw", "log", "cpm", "logcpm"),
  cols = NULL,
  pt.shape = 16,
  pt.size = 0.3,
  pt.alpha = 0.5,
  label = NULL,
  cols.scale = NULL,
  reverseY = NULL,
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- group.by:

  values to group points by. Must be in colData of spe. If NULL, will
  try with 'cols' if available.

- feature:

  Feature to group polygons by. Must be in rownames(spe).

- assay:

  Name of assay to use for plotting feature.

- type:

  Transformation to apply for the group/feature. Options are "raw" ,
  "log", "cpm", "logcpm", or a function that accepts and returns a
  vector of the same length.

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

- cols.scale:

  vector of position for color if colors should not be evenly
  positioned. See
  [scale_color_gradientn](https://ggplot2.tidyverse.org/reference/scale_gradient.html).
  Only applicable for continuous values.

- reverseY:

  Logical. Whether to reverse Y coordinates. Default is TRUE if the spe
  contains an image (even if not plotted) and FALSE if otherwise.

- ...:

  Parameters pass to plotImage

## Value

A ggplot object.

## Examples

``` r

data("xenium_bc_spe")

plotSpatial(spe, group.by = "cell_type", pt.size = 0.5, pt.alpha = 0.6)

```
