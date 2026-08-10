# Read Visium output into spe

Read Visium output into spe

## Usage

``` r
readVisium(
  dir,
  sample_id = "sample01",
  count = NULL,
  coord = NULL,
  image = NULL,
  scale_factors = NULL,
  feature_type = "Gene Expression",
  pixel_to_micron = TRUE
)
```

## Arguments

- dir:

  directory containing the Visium files

- sample_id:

  Name of the sample.

- count:

  Name of the h5 file with the count assay.

- coord:

  Path to the tissue coordinates file (csv or parquet), or a data.frame
  of coordinates (rownames = barcodes) with columns 'pxl_col_in_fullres'
  and 'pxl_row_in_fullres'.

- image:

  Names of the image files.

- scale_factors:

  Names of the scale factors file

- feature_type:

  Feature type to retain. Defaults to "Gene Expression" to exclude
  non-gene features. Set to NULL to keep all features.

- pixel_to_micron:

  Logical. If TRUE (default), convert the spot coordinates from
  full-resolution pixels to microns using 'spot_diameter_fullres' from
  the scale factors file, and store the conversion factor in
  metadata(spe)\$um_per_pixel. Set to FALSE to keep the coordinates in
  pixels (previous behaviour).
