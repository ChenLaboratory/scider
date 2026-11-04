# Plot model statistics using heatmap.

Plot model statistics using heatmap.

## Usage

``` r
plotCorHeatmap(
  model.result,
  stats = c("cor.coef", "cor.expected", "shift", "t", "p.Pos", "p.Neg"),
  roi = "all",
  cell.type = "all",
  silent = FALSE,
  level = c("overall", "ROI"),
  ...
)
```

## Arguments

- model.result:

  The output of corDensity(), a list with components 'ROI' and
  'overall', or one of these DataFrames.

- stats:

  Character value. The statistic to plot, one of "cor.coef" (default),
  "cor.expected", "shift", "t", "p.Pos" or "p.Neg". "cor.expected" and
  "shift" are only reported by corDensity(method = "permutation") and
  "t" only by corDensity(method = "modified.ttest").

- roi:

  Character value. By default is all. The specific ROIs to be plotted.

- cell.type:

  Character value. By default is all. The cell types to be plotted.

- silent:

  Do not draw the plot (useful when using the gtable output).

- level:

  Character value, either "overall" or "ROI". Which component to plot
  when 'model.result' is the list returned by corDensity(). Defaults to
  "overall", or to "ROI" when specific ROIs are selected with 'roi'.
  Ignored when 'model.result' is a DataFrame.

- ...:

  Other parameters passed to pheatmap().

## Value

A pheatmap object.

## Examples

``` r

data("xenium_bc_spe")

coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")

spe <- gridDensity(spe, coi = coi)

spe <- findROI(spe, coi = coi)

model_result <- corDensity(spe, roi = coi, nperm = 199)

plotCorHeatmap(model_result)


plotCorHeatmap(model_result, level = "ROI", stats = "shift")

```
