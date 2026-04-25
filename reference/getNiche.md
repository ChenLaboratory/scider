# Build a niche assay based on the profile of neighbouring cells

Build a niche assay based on the profile of neighbouring cells

## Usage

``` r
getNiche(
  spe,
  at = c("cell", "grid"),
  nbrs_name = NULL,
  group.by,
  use_weight = FALSE
)
```

## Arguments

- spe:

  A SpatialExperiment object

- at:

  Option of cell or grid neighbourhood

- nbrs_name:

  Name of the neighbour list in `spe@metadata$grid[[at]]`

- group.by:

  Character vector to group neighbours cell by. Should be in either
  colData(spe) or spe@metadata\$grid_density, depending on "at".
  Multiple groups can be used. See details

- use_weight:

  Whether to scale each nbr based on its weight

## Value

A matrix where rows are cells/grid points and cols are groups based on
group.by

## Details

For numerical group, result will be sum of nbrs for each cell. For
categorical group (factor/string), result will be counts of nbrs
belonging in category

## Examples

``` r
data("xenium_bc_spe")

spe <- findNbrsSpatial(spe, k=30)
niche <- getNiche(spe, at="cell", group.by="cell_type")
```
