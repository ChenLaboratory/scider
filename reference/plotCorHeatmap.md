# Plot model statistics using heatmap.

Plot model statistics using heatmap.

## Usage

``` r
plotCorHeatmap(
  model.result,
  stats = c("cor.coef", "t", "p.Pos", "p.Neg"),
  roi = "all",
  cell.type = "all",
  silent = FALSE
)
```

## Arguments

- model.result:

  A dataFrame object.

- stats:

  Character value. Choose either coefficient or t. Coefficient by
  default.

- roi:

  Character value. By default is all. The specific ROIs to be plotted.

- cell.type:

  Character value. By default is all. The cell types to be plotted.

- silent:

  Do not draw the plot (useful when using the gtable output).

## Value

A pheatmap object.

## Examples

``` r

data("xenium_bc_spe")

coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")

spe <- gridDensity(spe, coi = coi)

spe <- findROI(spe, coi = coi)

model_result <- corDensity(spe, roi = coi)

plotCorHeatmap(model_result$ROI)

```
