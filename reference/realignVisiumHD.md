# Scale and straighten out VisiumHD coordinates

Scale and straighten out VisiumHD coordinates

## Usage

``` r
realignVisiumHD(spe, distPoint = NULL)
```

## Arguments

- spe:

  A SpatialExperiment object.

- distPoint:

  Numeric. Desired point to point distance. If NULL, will try to
  determine the bin level of spe and use that.

## Details

This function rescale the distance between points to 8um (or other
value) to matchthe real distance. In addition, Visium spots have a
slight tilt to them which this function will also fix
