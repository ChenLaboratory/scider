# UMAP using uwot. Parameters are set to be similar to Seurat's

UMAP using uwot. Parameters are set to be similar to Seurat's

## Usage

``` r
runUMAP(
  spe,
  n_neighbors = 30,
  n_components = 2,
  metric = "cosine",
  min_dist = 0.3,
  assay = NULL,
  dimred = "PCA",
  n_dimred = NULL,
  name = "UMAP",
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- n_neighbors, n_components, metric, min_dist:

  See [umap](https://jlmelville.github.io/uwot/reference/umap.html)

- assay:

  Name of assay for UMAP. Incompatible with dimred.

- dimred:

  Name of the dimensionality reduction (e.g. PCA) for UMAP. Incompatible
  with assay

- n_dimred:

  Integer scalar or vector specifying the dimensions to use if dimred is
  specified.

- name:

  Name to store the UMAP in the spe's
  [reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html).

- ...:

  Other parameters to be passed to
  [umap](https://jlmelville.github.io/uwot/reference/umap.html).

## Value

A SpatialExperiment with the UMAP stored in
[reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html).

## Details

By default, runUMAP uses PCA (from
[runPCA](https://chenlaboratory.github.io/scider/reference/runPCA.md)).
If that's unavailable, it falls back to logcounts, then counts assay.

## Examples

``` r

data("xenium_bc_spe")
spe <- runPCA(spe)
#> Default assay logcounts not found. Switching to counts assay instead.
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
spe <- runUMAP(spe,dimred="PCA",n_dimred=10)
```
