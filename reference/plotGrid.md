# Plot grid from metadata.

Plot grid from metadata.

## Usage

``` r
plotGrid(
  spe,
  group.by = NULL,
  feature = NULL,
  assay = "counts",
  type = c("raw", "log", "cpm", "logcpm"),
  cols = NULL,
  pol.border = FALSE,
  pol.alpha = 1,
  probs = 0,
  cutoff = NULL,
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

  values to group polygons by. Must be in spe@metadata\$grid_density, or
  colData(spe) if gridLevelAnalysis is TRUE. If NULL, will try with cols
  if available.

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

- pol.border:

  Boolean. Whether to draw border for each polygon.

- pol.alpha:

  alpha of points between 0 and 1.

- probs:

  Numeric value between 0 and 1, used for filtering uninformative grid.
  Only applicable for continuous values.

- cutoff:

  Numeric. Either a vector of length 2 for the lower & upper bounds of
  data to be included, or length 1 for the lower bound. Override probs
  if specified. Only applicable for continuous values.

- label:

  label for the legend

- cols.scale:

  vector of position for color if colors should not be evenly
  positioned. See
  [scale_fill_gradientn](https://ggplot2.tidyverse.org/reference/scale_gradient.html).
  Only applicable for continuous values.

- reverseY:

  Logical. Whether to reverse Y coordinates. Default is TRUE if the spe
  contains an image (even if not plotted) and FALSE if otherwise.

- ...:

  Parameters pass to
  [plotImage](https://chenlaboratory.github.io/scider/reference/plotImage.md)

## Value

A ggplot object.

## Examples

``` r

data("xenium_bc_spe")

spe <- gridDensity(spe)

plotGrid(spe, group.by = "density_overall")

```
