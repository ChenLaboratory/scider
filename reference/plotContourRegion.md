# Visualising an sf object (for internal use only at the moment)

Visualising an sf object (for internal use only at the moment)

## Usage

``` r
plotContourRegion(
  spe,
  coi,
  id = "cell_type",
  overlay = c("density", "cell"),
  sub.level
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- coi:

  A character vector of length 1 of the cell type of interest.

- id:

  A character. The name of the column of colData(spe) containing the
  cell type identifiers. Set to cell_type by default.

- overlay:

  Character vector. Either plot overlay on density or cells.

- sub.level:

  Numeric vector of length 1 or 2, identifies which density level to
  plot. When length is 1, plot the density region above this level. When
  length is 2, plot the density region between the two levels.

## Value

A ggplot object.
