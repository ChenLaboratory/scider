# Check which cells are in which regions

Check which cells are in which regions

## Usage

``` r
cellsInRegion(spe, region, name_to, NA_level = "0", levels = NULL)
```

## Arguments

- spe:

  A SpatialExperiment object.

- region:

  List or an sf object that represents a region or an ROI.

- name_to:

  Colname in colData(spe) to store the annotation.

- NA_level:

  Label for cells not falling in any of the regions. Default to 0.

- levels:

  Factor levels.

## Value

A SpatialExperiment object. The region information of each cell is
stored in the colData.
