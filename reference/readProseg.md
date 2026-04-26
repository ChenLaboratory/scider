# Read Proseg V2 output into spe

Read Proseg V2 output into spe

## Usage

``` r
readProseg(
  dir,
  sample_id = "sample01",
  count = "expected-counts\\.(csv|mtx)\\.gz",
  coord = "cell-metadata.csv.gz",
  gene = "gene-metadata.csv.gz",
  coordNames = c("centroid_x", "centroid_y", "centroid_z")
)
```

## Arguments

- dir:

  directory containing the Proseg files

- sample_id:

  Name of the sample.

- count:

  Name of the file with the count assay.

- coord:

  Name of the file with the tissue coordinates

- gene:

  Name of the file with the gene metadata. This

- coordNames:

  Name of the coordinates for the spe

## Details

This does not work on zarr file output of proseg V3
