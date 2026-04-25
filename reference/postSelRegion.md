# Merge sel_region from the selectRegion function to SpatialExperiment.

Merge sel_region from the selectRegion function to SpatialExperiment.

## Usage

``` r
postSelRegion(spe, sel_region)
```

## Arguments

- spe:

  A SpatialExperiment object.

- sel_region:

  A dataframe object. Can be generated from function selectRegion.

## Value

A SpatialExperiment object.

## Examples

``` r
data("xenium_bc_spe")

coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")

spe <- gridDensity(spe, coi = coi)

sel_region <- data.frame(
    "node" = seq(10),
    "node_x" = seq(10),
    "node_y" = seq(10)
)

spe1 <- postSelRegion(spe, sel_region)
```
