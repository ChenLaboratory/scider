# Combine the top markers of each cluster into one data frame.

Convenience wrapper around
[getMarkers](https://chenlaboratory.github.io/scider/reference/getMarkers.md)
output: takes the top `n` genes from each cluster's (already ranked)
table and stacks them into a single data frame, with a leading `cluster`
column and a `gene` column (the same gene can be a top marker for more
than one cluster, so row names cannot stay unique). Analogous to the
single table returned by `Seurat::FindAllMarkers`.

## Usage

``` r
topMarkers(markers, n = 10, p.value = 1, min.prop = 0)
```

## Arguments

- markers:

  The named list returned by
  [getMarkers](https://chenlaboratory.github.io/scider/reference/getMarkers.md).

- n:

  Integer. Number of top genes to take from each cluster. Default 10.

- p.value:

  Numeric. Keep only genes whose adjusted p-value (or raw `PValue` if no
  adjusted column is present) is at or below this cutoff, applied before
  taking the top `n`. Default 1 (no filtering).

- min.prop:

  Numeric in \[0, 1\]. Drop genes whose detection rate (the `Prop`
  column - proportion of cells expressing in the cluster) is below this
  value before taking the top `n`, so the slots are back-filled from
  further down the ranking. Default 0 (no filtering).

## Value

A data frame with columns `cluster`, `gene`, and the per-cluster
statistics from
[getMarkers](https://chenlaboratory.github.io/scider/reference/getMarkers.md),
ordered by cluster then by each cluster's ranking.

## Examples

``` r

data("xenium_bc_spe")
spe <- normalizeAssay(spe)
spe <- runPCA(spe)
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
spe <- findNbrsSNN(spe, dimred = "PCA")
spe <- getClusters(spe, resolution = 0.5)
markers <- getMarkers(spe)
top <- topMarkers(markers, n = 10)
```
