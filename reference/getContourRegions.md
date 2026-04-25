# Calculate areas between every two density levels

Calculate areas between every two density levels

## Usage

``` r
getContourRegions(spe, contour_name)
```

## Arguments

- spe:

  A SpatialExperiment object.

- contour_name:

  Name of contour in spe@metadata

## Value

A list of sf objects, each representing the region between two contour
density levels.
