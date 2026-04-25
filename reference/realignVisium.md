# Scale and straighten out Visium coordinates

Scale and straighten out Visium coordinates

## Usage

``` r
realignVisium(spe, distPoint = 100)
```

## Arguments

- spe:

  A SpatialExperiment object.

- distPoint:

  Numeric. Desired point to point distance.

## Details

This function rescale the distance between points to 100um (or other
value) to matchthe real distance. In addition, Visium spots have a
slight tilt to them which this function will also fix
