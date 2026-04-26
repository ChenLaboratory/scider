# Read Xenium output into spe

Read Xenium output into spe

## Usage

``` r
readXenium(
  dir,
  sample_id = "sample01",
  count = file.path(dir, "cell_feature_matrix.h5"),
  coord = file.path(dir, "cells.parquet"),
  image = file.path(dir, "morphology.ome.tif"),
  image_reso = 6,
  image_layer = NULL
)
```

## Arguments

- dir:

  directory containing the Xenium files

- sample_id:

  Name of the sample.

- count:

  Name of the h5 file with the count assay.

- coord:

  Name of the parquet file with the tissue coordinates

- image:

  Names of the ome.tif image files.

- image_reso:

  resolution of the image to use. From 1-8 (lower = better resolution).
  See https://kb.10xgenomics.com/hc/articles/11636252598925. Default to
  6

- image_layer:

  Which layer of the tiff image to use. Default is the middle-most layer
