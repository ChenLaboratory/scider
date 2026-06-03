# Test for density correlation between two cell types.

Test for density correlation between two cell types.

## Usage

``` r
corDensity(spe, coi = NULL, roi = NULL, probs = 0.85, trace = FALSE)
```

## Arguments

- spe:

  A SpatialExperiment object.

- coi:

  Character vector for cell types of interest for density correlation
  analysis. Default is NULL, which is to consider all cell types
  previously calculated in the gridDensity() step.

- roi:

  Character. The name of the group or cell type on which the roi is
  computed. Default is NULL for no subsetting cell types by ROI

- probs:

  A numeric scalar. The threshold of proportion that used to filter
  grids by density when ROIs have not been identified previously.
  Ignored if 'roi' is present in the 'metadata' component of spe.
  Default to 0.85.

- trace:

  Logical. If TRUE, print the process of testing. Default to FALSE.

## Value

A DataFrame containing the testing results.

## Examples

``` r

data("xenium_bc_spe")

coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")

spe <- gridDensity(spe, coi = coi)

spe <- findROI(spe, coi = coi, method = "walktrap")

result <- corDensity(spe, roi = coi)
```
