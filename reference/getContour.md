# Get contour from density

Get contour from density

## Usage

``` r
getContour(
  spe,
  coi = NULL,
  equal.cell = TRUE,
  bins = NULL,
  binwidth = NULL,
  breaks = NULL,
  id = NULL
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- coi:

  A character vector of cell types of interest (COIs). All cell types
  are chosen if NULL or `overall`.

- equal.cell:

  Logical. Whether to use produce contour levels so that there are
  roughly the same number of cells of the COI at each level. Default to
  TRUE.

- bins:

  An integer. Number of contour levels.

- binwidth:

  A numeric scale of the smoothing bandwidth.

- breaks:

  A numeric scale referring to the breaks in `ggplot2:::contour_breaks`.

- id:

  A character. The name of the column of colData(spe) containing the
  cell type identifiers. Set to cell_type by default or in_tissue if spe
  is Visium. Only needed when `equal.cell = TRUE`.

## Value

A SpatialExperiment object. An sf object of the contour region of the
specified level is stored in the metadata of the SpatialExperiment
object.

## Examples

``` r

data("xenium_bc_spe")

spe <- gridDensity(spe)

coi <- "Breast cancer"

spe <- getContour(spe, coi = coi)
#> Using bins = 10 to draw contours with equal cell numbers.
```
