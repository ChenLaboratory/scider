# Annotate all cells with contour level of cell type-specific density.

Annotate all cells with contour level of cell type-specific density.

## Usage

``` r
allocateCells(
  spe,
  to.roi = TRUE,
  roi = NULL,
  to.contour = TRUE,
  contour = NULL
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- to.roi:

  Logical. Whether to allocate cells to ROIs.

- roi:

  Character. The name of the group or cell type on which the roi is
  computed. If NULL, then the cell allocation will be performed for all
  detected roi Default to NULL.

- to.contour:

  Logical. Whether to allocate cells to contour levels.

- contour:

  Character. The name of the group or cell type on which the contour
  level is computed. If NULL, then the cell allocation will be performed
  for all detected contours. Default to NULL.

## Value

A SpatialExperiment object. An extra column is added to the colData.

## Examples

``` r
data("xenium_bc_spe")
spe <- gridDensity(spe)
coi <- "Breast cancer"
spe <- findROI(spe, coi = coi)
spe <- getContour(spe, coi = coi)
#> Using bins = 10 to draw contours with equal cell numbers.
spe <- allocateCells(spe, contour = coi)
#> Assigning cells to ROIs defined by Breast cancer 
#> Assigning cells to contour levels of Breast cancer 
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
```
