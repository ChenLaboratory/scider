# Renumber all clusters by size.

Renumbers the clusters in a colData column `1..K` by size (largest
first), with "unassigned" kept as the last level. Useful for tidying up
the compound labels left by `getSubClusters(relabel = FALSE)` /
`mergeClusters(relabel = FALSE)` after a multi-step workflow.

## Usage

``` r
relabelClusters(
  spe,
  cluster_name = "cluster",
  new_name = cluster_name,
  start_from = NULL
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- cluster_name:

  Name of the cluster column in
  [colData](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  to renumber.

- new_name:

  Name of the column to store the result. Defaults to overwriting
  cluster_name.

- start_from:

  Integer at which numbering starts. Defaults to NULL, which infers the
  base (0- or 1-based) from the existing labels.

## Value

A SpatialExperiment with the renumbered clusters in colData.

## Details

Renumbering is size-based, so run it **before** cell-type annotation -
relabelling after a cluster-to-cell-type mapping would scramble that
mapping.

## Examples

``` r

data("xenium_bc_spe")
spe <- normalizeAssay(spe)
spe <- runPCA(spe)
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
spe <- findNbrsSNN(spe, dimred = "PCA")
spe <- getClusters(spe, resolution = 0.5)
spe <- getSubClusters(spe, cluster = 2, resolution = 0.5)
# Tidy the 2_1/2_2/... labels into plain 1..K numbers by size
spe <- relabelClusters(spe)
```
