# Construct a SNN neighbour list from assay.

Construct a SNN neighbour list from assay.

## Usage

``` r
findNbrsSNN(
  spe,
  assay = NULL,
  dimred = "PCA",
  n_dimred = 10,
  k = 20,
  BNPARAM = BiocNeighbors::AnnoyParam(),
  type = c("rank", "number", "jaccard"),
  nbrs_name = NULL,
  cpu_threads = 6
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- assay:

  Name of assay for clustering. Incompatible with dimred.

- dimred:

  Name of the dimensionality reduction (e.g. PCA) for clustering.
  Incompatible with assay

- n_dimred:

  Integer scalar or vector specifying the dimensions to use if dimred is
  specified.

- k:

  Integer scalar for number of nearest neighbors to find.

- BNPARAM:

  [BiocNeighborParam](https://rdrr.io/pkg/BiocNeighbors/man/BiocNeighborParam.html)
  object specifying the nearest neighbor algorithm. Default is Annoy.

- type:

  Type of weighting scheme for shared neighbors. Options are rank,
  number, and jaccard. type="rank" is defined in Xu and Su (2015).

- nbrs_name:

  Name of the neighbour list to be stored in spe. Default to be
  assay/dimred + "\_snn".

- cpu_threads:

  Number of cpu threads for parallel computation.

## Value

A spe with the clusters stored in
[reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html).

A SpatialExperiment object

## Details

Construct a SNN neighbour list using either the spe's assay or reduced
dimension and store it in spe@metadata\$nbrs\$cell

neighbour list contain

## Examples

``` r
data("xenium_bc_spe")
spe <- runPCA(spe)
#> Default assay logcounts not found. Switching to counts assay instead.
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
spe <- findNbrsSNN(spe,dimred="PCA")
#> [1] "Getting K-nearest neighbour"
#> [1] "Getting shared nearest neighbour"
```
