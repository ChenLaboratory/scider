# Fast PCA using irlba.

Fast PCA using irlba.

## Usage

``` r
runPCA(
  spe,
  n_pcs = 50,
  assay = "logcounts",
  centre = TRUE,
  scale = TRUE,
  clip = NULL,
  name = "PCA",
  genes = "hvg",
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- n_pcs:

  Number of principal components to calculate

- assay:

  Name of assay used for PCA. See details for defaults.

- centre:

  Logical. Whether to centre the assay before PCA.

- scale:

  Logical. Whether to scale the variance to 1 before PCA.

- clip:

  Maximum absolute z-score after scaling. Values beyond this are
  clipped. Prevents rare-marker genes (near-zero SD) from creating
  extreme outliers that fragment the UMAP. Defaults to `10`. Set to
  `Inf` to disable.

- name:

  Name to store the PCA in the spe's
  [reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)

- genes:

  Subset of features for PCA. Can be a column in rowData or a vector of
  gene names, indices, or booleans. Default to hvg if
  [getHVG](https://chenlaboratory.github.io/scider/reference/getHVG.md)
  was run.

- ...:

  Other parameters to be passed to
  [irlba](https://rdrr.io/pkg/irlba/man/irlba.html).

## Value

A SpatialExperiment with the PCA stored in
[reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html).

## Details

By default, runPCA uses logcounts assay (from
[normalizeAssay](https://chenlaboratory.github.io/scider/reference/normalizeAssay.md)).
If that's unavailable, it falls back to counts assay

## Examples

``` r

data("xenium_bc_spe")

spe <- runPCA(spe)
#> Default assay logcounts not found. Switching to counts assay instead.
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
```
