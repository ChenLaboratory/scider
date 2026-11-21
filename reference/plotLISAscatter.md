# Scatterplot for local moran's I

plot result obtained from localMoran()

## Usage

``` r
plotLISAscatter(
  lisa,
  quadrant.count = TRUE,
  text.size = 11/.pt,
  xlab = "Data",
  ylab = "Spatial Lag"
)
```

## Arguments

- lisa:

  A list obtained from
  [localMoran](https://chenlaboratory.github.io/scider/reference/localMoran.md)

- quadrant.count:

  Whether to count values at each quadrant (Low-Low, Low-High,
  High-High, High-Low)

- text.size:

  Numeric for text size of quadrant.count

- xlab:

  label for the x-axis

- ylab:

  label for the y-axis

## Value

A ggplot object.

## Details

Each point is a cell or grid point, with its standardised value on the
x-axis and its spatial lag on the y-axis. The title gives the global
Moran's I, as returned by
[globalMoran](https://chenlaboratory.github.io/scider/reference/globalMoran.md).
The line is the least-squares fit of the spatial lag on the value; its
slope equals Moran's I only when the weights are row-standardised.

## Examples

``` r

data("xenium_bc_spe")
dat <- spe$total_counts
spe <- findNbrsSpatial(spe,k=50)
res <- localMoran(spe,data1=dat,at="cell")
plotLISAscatter(res)
```
