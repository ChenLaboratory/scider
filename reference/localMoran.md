# Calculate local Moran for 1 to 2 variables.

Calculate local Moran for 1 to 2 variables.

## Usage

``` r
localMoran(
  spe,
  data1,
  data2 = data1[],
  at = c("grid", "cell"),
  nbrs_name = NULL,
  hhonly = FALSE,
  significance_cutoff = 0.05,
  permutations = 999,
  seed = 123456789,
  cpu_threads = 6
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- data1:

  Numeric vector 1. Must be same length as number of grid points or
  cells, depending on 'at'.

- data2:

  Numeric vector 2 for bivariate local Moran. Must be same length as
  data1.

- at:

  Option of grid or cell for where to look for neighbour list

- nbrs_name:

  Name of the neighbour list in `spe@metadata$nbrs[[at]]` for Moran's I.
  If NULL, uses the most recent spatial list from
  [findNbrsSpatial](https://chenlaboratory.github.io/scider/reference/findNbrsSpatial.md)
  at cell level, or the most recent grid list from
  [findNbrsGrid](https://chenlaboratory.github.io/scider/reference/findNbrsGrid.md)
  at grid level.

- hhonly:

  Only high-high clusters, which is more interpretable. Other clusters
  (e.g. "Low-Low","High-Low",...) will be assigned as undefined.

- significance_cutoff:

  Cutoff for p-value to filter non-significant clusters

- permutations:

  Number of permutations for p-value.

- seed:

  Integer. For random permutations.

- cpu_threads:

  The number of cpu threads used for parallel computation.

## Value

List with the local Moran values (lisa), clusters (cluster), pseudo
p-values (p.value), spatial lags (spatiallag) and the sum of all
neighbour weights (W). `sum(lisa)/W` is the global Moran's I, as
returned by
[globalMoran](https://chenlaboratory.github.io/scider/reference/globalMoran.md).

## Examples

``` r
data("xenium_bc_spe")

## At grid.
spe <- gridDensity(spe, coi = "Breast cancer")
dat <- spe@metadata$grid_density$density_breast_cancer
spe <- findNbrsGrid(spe)
res <- localMoran(spe,data1=dat, at = "grid")

## At cell.
dat <- as.numeric(spe$cell_type=="Breast cancer")
spe <- findNbrsSpatial(spe,k=20)
res <- localMoran(spe,data1=dat,at="cell")
```
