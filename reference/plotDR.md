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
  type = c("log", "raw", "cpm", "logcpm"),
  cols = NULL,
  highlight = NULL,
  cols.highlight = NULL,
  pt.shape = 16,
  pt.size = 0.7,
  pt.size.highlight = 1,
  pt.alpha = 0.6,
  label = NULL,
  xlab = NULL,
  ylab = NULL,
  cols.scale = NULL,
  ncol = NULL,
  per.scale = TRUE,
  range = NULL,
  transform = "identity"
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

  Feature(s) to colour points by; must be in rownames(spe). If a vector
  of more than one feature is supplied, one reduced-dimension panel is
  drawn per feature, faceted (see `ncol`), with a shared colour scale.

- assay:

  Name of assay to use for plotting feature. Default "counts".

- type:

  Transformation applied to feature expression: "log" (default,
  log2(1+x)), "raw" (no transform), "cpm", or "logcpm" (log2 CPM). For
  multiple features the legend title defaults to the matching unit
  ("log2 Cts", "Counts", "CPM", "log2-CPM"); for a single feature the
  unit is the legend title and the gene name becomes the plot title.
  Override the legend with `label`.

- cols:

  Colour palette. Can be a vector of colours or a function that accepts
  an integer n and return n colours.

- highlight:

  Optional cells to emphasise, given as either a vector of group.by
  levels (characters or cluster numbers), or a logical vector of length
  ncol(spe) selecting cells directly. Highlighted cells are drawn last
  (on top) at pt.size.highlight; all other cells are light grey at
  pt.size.

- cols.highlight:

  Colour(s) for the highlighted cells. Defaults to NULL, which keeps
  each level's usual group.by palette colour. A single colour (e.g.
  "red") colours all highlighted cells the same; a vector matching the
  number of 'highlight' entries gives one colour per level (matched by
  position).

- pt.shape:

  shape of points.

- pt.size:

  size of points.

- pt.size.highlight:

  size of highlighted points (see highlight).

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

- ncol:

  Number of columns when plotting multiple features. Passed to
  [facet_wrap](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  (shared scale) or
  [wrap_plots](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
  (per-panel scales). Default NULL lets the layout be chosen
  automatically.

- per.scale:

  Logical. For multiple features, whether each panel gets its own colour
  scale (TRUE, default; like `Seurat::FeaturePlot`, via the 'patchwork'
  package) or a single shared colour scale across panels (FALSE).

- range:

  Numeric length-2 (lower, upper) cap for continuous colour values

  - a single feature, or multiple features with a shared scale; values
    outside are clamped. A single value is taken as the upper bound.
    Default NULL (no capping). Ignored when `per.scale = TRUE` (each
    panel auto-scales).

- transform:

  Name of a transformation for the continuous colour scale (e.g.
  "log10", "log1p", "pseudo_log", "sqrt"), passed to
  [scale_color_gradientn](https://ggplot2.tidyverse.org/reference/scale_gradient.html);
  the colour spectrum is spaced by the transform while the legend stays
  in original units. Default "identity" (no transformation). Use "log10"
  for nicely log-spaced legend breaks (needs positive values);
  "log1p"/"pseudo_log" tolerate zeros but keep linear breaks. Only
  affects continuous (feature or numeric group.by) colouring.

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
