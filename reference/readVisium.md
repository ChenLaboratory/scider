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
  scale_factors = NULL
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

  Name of the csv file with the tissue coordinates

- image:

  Names of the image files.

- scale_factors:

  Names of the scale factors file
