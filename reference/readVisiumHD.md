# Read VisiumHD output into spe

Read VisiumHD output into spe

## Usage

``` r
readVisiumHD(dir, bin = c("016um", "008um", "002um", "segmented"), ...)
```

## Arguments

- dir:

  directory containing the VisiumHD files

- bin:

  Which output to read. Bin sizes "016um", "008um", "002um" read the
  corresponding 'binned_outputs/square\_\*' folder. "segmented" reads
  the cell-segmentation results in 'segmented_outputs': the count matrix
  'filtered_feature_cell_matrix.h5', with per-cell coordinates taken
  from the centroids of 'cell_segmentations.geojson'.

- ...:

  Parameters for readVisium

## Details

For "segmented", cells have no array_row/array_col grid, so the result
is cell-level (like Xenium) and is not compatible with the Visium
spot-grid options of gridDensity() / trimEdge().
