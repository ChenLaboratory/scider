# Summarize a SpatialExperiment object at grid-level

Summarize a SpatialExperiment object at grid-level

## Usage

``` r
gridSPE(spe, cell.count = FALSE, id = "cell_type", split.count.by = id)
```

## Arguments

- spe:

  A SpatialExperiment object.

- cell.count:

  Logical. Whether to obtain the number of cells within each group
  identified by the 'id' column in colData(spe). Default to FALSE.

- id:

  A character. The name of the column of colData(spe) containing the
  cell type identifiers. Set to 'cell_type' by default.

- split.count.by:

  A character. The name of the column of colData(spe). When it is not
  NULL, a grid-level count matrix is calculated for each member
  specified in that column of colData(spe) and stored in the
  assays(spe). Set to 'cell_type' by default.

## Value

A SpatialExperiment object.

## Examples

``` r

data("xenium_bc_spe")

spe <- gridDensity(spe)

spe_grid <- gridSPE(spe)
```
