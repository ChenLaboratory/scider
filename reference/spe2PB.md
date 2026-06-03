# Given a 'SpatialExperiment' data object, create pseudo-bulk samples using the colData information and return a DGEList object

Given a 'SpatialExperiment' data object, create pseudo-bulk samples
using the colData information and return a DGEList object

## Usage

``` r
spe2PB(
  spe,
  by.group = TRUE,
  group.id = "cell_type",
  keep.groups = NULL,
  roi = NULL,
  roi.only = TRUE,
  contour = NULL
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- by.group:

  Logical. Whether to perform pseudo-bulking by group. TRUE by default.

- group.id:

  Character. The column name of the colData(spe) that contains the group
  information. Default to 'cell_type'.

- keep.groups:

  Vector. Values from group.id to include in pseudo- bulking. Default is
  NULL, where all cells are included in pseudo-bulking.

- roi:

  Character. The name of the group or cell type on which the roi is
  computed. If NULL, then no pseudo-bulking will be performed based on
  roi. Default to NULL.

- roi.only:

  Logical. Whether to filter out pseudo-bulk samples formed by cells not
  in any ROIs. TRUE by default.

- contour:

  Character. The name of the group or cell type on which the contour
  level is computed. If NULL, then no pseudo-bulking will be performed
  based on contour level. Default to NULL.

## Value

An edgeR::DGEList object where each library (column) is a pseudo-bulk
sample.

## Examples

``` r

data("xenium_bc_spe")

spe <- gridDensity(spe)

coi <- "Breast cancer"

spe <- findROI(spe, coi = coi)

spe <- allocateCells(spe, to.contour=FALSE)
#> Assigning cells to ROIs defined by Breast cancer 

y <- spe2PB(spe, roi = coi)
```
