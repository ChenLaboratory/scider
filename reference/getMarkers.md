# Find up-regulated marker genes for all clusters (quasi-NB GLM).

Pseudo-bulks cells by cluster via
[spe2PB](https://chenlaboratory.github.io/scider/reference/spe2PB.md)
and, for each cluster, tests its genes against the average of all other
clusters (1-vs-rest) using a quasi-likelihood approach. Because the full
design is saturated (one pseudo-bulk per cluster, 0 residual df), each
contrast is tested by fitting the corresponding null design (one fewer
column, 1 residual df) with
[glmQLFit](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) and using the
quasi-deviance chi-square statistic. The test is **one-sided**: p-values
reflect up-regulation in the cluster, so each table ranks that cluster's
positive markers.

## Usage

``` r
getMarkers(
  spe,
  cluster_name = "cluster",
  dispersion = 0.05,
  method = c("quasi", "lrt"),
  adjust.method = "BH",
  sort.by = "PValue",
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- cluster_name:

  Character. Column name in `colData(spe)` containing cluster labels.
  Default `"cluster"`.

- dispersion:

  Numeric. Fixed NB dispersion (BCV^2) used for the GLM fits. Required
  because the saturated full design (one pseudo-bulk per cluster) leaves
  no residual df for dispersion estimation. Default 0.05 (BCV
  \\\approx\\ 0.224).

- method:

  Character. Testing pipeline: "quasi" (default) fits each null design
  with [glmQLFit](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) and uses
  the quasi-deviance chi-square statistic; "lrt" uses the standard
  [glmFit](https://rdrr.io/pkg/edgeR/man/glmfit.html) +
  [glmLRT](https://rdrr.io/pkg/edgeR/man/glmLRT.html) likelihood-ratio
  test. Both use the preset `dispersion`.

- adjust.method:

  Character. Multiple-testing adjustment method passed to
  [p.adjust](https://rdrr.io/r/stats/p.adjust.html) (e.g. "BH", "BY",
  "holm", "bonferroni", "none"). Default "BH".

- sort.by:

  Character. How to order each table: "PValue" (default), "logFC" (by
  absolute log fold change), "stat", or "none".

- ...:

  Additional arguments passed to
  [glmQLFit](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) for the
  null-model fits (used only when `method = "quasi"`).

## Value

A named `list` of data frames, one per cluster, each containing all
genes. Columns:

- logFC:

  log2 fold change of the cluster vs the average of the rest.

- logCPM:

  average log2 counts-per-million across pseudo-bulk samples.

- Prop:

  proportion of cells in the cluster with a non-zero count for the gene.

- stat:

  signed quasi statistic (signed square root of the quasi chi-square
  statistic); positive means up-regulated in the cluster.

- PValue:

  one-sided (up-regulation) p-value.

- FDR:

  adjusted p-value (column named per `adjust.method`; "FDR" for BH/BY,
  "FWER" for family-wise methods). Absent when `adjust.method = "none"`.

## Details

Each per-cluster table holds **all** genes, sorted by `sort.by` with an
adjusted p-value column added via `adjust.method`. Use
[topMarkers](https://chenlaboratory.github.io/scider/reference/topMarkers.md)
to extract and combine the top genes per cluster. For a direct
comparison between two clusters (or two groups of clusters), see
[getDE](https://chenlaboratory.github.io/scider/reference/getDE.md).

## Examples

``` r

data("xenium_bc_spe")
spe <- normalizeAssay(spe)
spe <- runPCA(spe)
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
spe <- findNbrsSNN(spe, dimred = "PCA")
#> [1] "Getting K-nearest neighbour"
#> [1] "Getting shared nearest neighbour"
spe <- getClusters(spe, resolution = 0.5)
#> Reassigned all 2 cells from clusters with <= 1 cells (0 by graph, 2 by distance).
markers <- getMarkers(spe)
```
