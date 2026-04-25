# Manually merge ROIs

Manually merge ROIs

## Usage

``` r
mergeROI(
  spe,
  roi = NULL,
  merge.list = NULL,
  remove.ids = NULL,
  id = "component",
  rename = FALSE
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- roi:

  Character. The name of the group or cell type on which the roi is
  computed. All cell types are chosen if NULL or 'overall'.

- merge.list:

  A (named) list of vectors of ROI ids to be merged. Each vector in the
  list should be of length greater than or equal to 2. If no name is
  specified, the merged ROI will be named by concatenating ROIs being
  merged.

- remove.ids:

  Optional. A vector of ROI ids to be removed.

- id:

  Character. The name of the column in `spe@metadata$roi` that stores
  the ROIs to be merged. Default is "component".

- rename:

  Logical. If TRUE, names of merge.list are ignored. ROIs will be given
  a new name. For the unmerged ROIs, their new names are not necessarily
  the same as those before merging.

## Value

A SpatialExperiment object.

## Examples

``` r
data("xenium_bc_spe")

coi <- c("Breast cancer", "Fibroblasts")

spe <- gridDensity(spe, coi = coi)

spe <- findROI(spe, coi = coi, method = "walktrap")

spe <- mergeROI(spe, roi = coi, list("1-2" = 1:2))
```
