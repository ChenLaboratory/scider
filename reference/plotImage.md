# Plot background image of spe

Plot background image of spe

## Usage

``` r
plotImage(
  spe,
  image = TRUE,
  image_id = NULL,
  sample_id = NULL,
  reverseY = NULL,
  crop = TRUE,
  image.alpha = 1
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- image:

  Logical. Whether to plot image (if present).

- image_id, sample_id:

  sample and image identifiers for scaling factor. See
  [scaleFactors](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html)
  for default behavior.

- reverseY:

  Logical. Whether to reverse Y coordinates. Default is TRUE if the spe
  contains an image (even if not plotted) and FALSE if otherwise.

- crop:

  Whether to crop the plot to the spots.

- image.alpha:

  alpha of points between 0 and 1.

## Value

a ggplot object if there is a valid image, else NULL.
