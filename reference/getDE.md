# Differential expression between two clusters or two groups of clusters.

Pseudo-bulks cells by cluster via
[spe2PB](https://chenlaboratory.github.io/scider/reference/spe2PB.md)
and performs a **two-sided** quasi-likelihood test comparing the average
of one set of clusters (`cluster`) against another (`cluster_2`),
analogous to
[glmQLFTest](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html). The
contrast is the difference of the two group means, tested by fitting the
null design with [glmQLFit](https://rdrr.io/pkg/edgeR/man/glmQLFit.html)
(see
[getMarkers](https://chenlaboratory.github.io/scider/reference/getMarkers.md)
for why the null design is used). The result table is post-processed
like [topTags](https://rdrr.io/pkg/edgeR/man/topTags.html).

## Usage

``` r
getDE(
  spe,
  cluster,
  cluster_2,
  cluster_name = "cluster",
  dispersion = 0.05,
  method = c("quasi", "lrt"),
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue",
  p.value = 1,
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- cluster:

  Character vector. One or more cluster labels forming the first group.
  Positive `logFC` means up-regulated in this group.

- cluster_2:

  Character vector. One or more cluster labels forming the second
  (reference) group. Must not overlap `cluster`.

- cluster_name:

  Character. Column name in `colData(spe)` containing cluster labels.
  Default `"cluster"`.

- dispersion:

  Numeric. Fixed NB dispersion (BCV^2) used for the GLM fits. Default
  0.05 (BCV \\\approx\\ 0.224).

- method:

  Character. Testing pipeline: "quasi" (default,
  [glmQLFit](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) quasi-deviance
  test) or "lrt" ([glmFit](https://rdrr.io/pkg/edgeR/man/glmfit.html) +
  [glmLRT](https://rdrr.io/pkg/edgeR/man/glmLRT.html)). Both use the
  preset `dispersion`.

- n:

  Integer. Maximum number of genes to return. Default `Inf` (all genes).

- adjust.method:

  Character. Multiple-testing adjustment method passed to
  [p.adjust](https://rdrr.io/r/stats/p.adjust.html). Default "BH".

- sort.by:

  Character. How to order the table: "PValue" (default), "logFC",
  "stat", or "none".

- p.value:

  Numeric. Cutoff on the adjusted p-value; genes above it are dropped.
  Default 1 (no filtering).

- ...:

  Additional arguments passed to
  [glmQLFit](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) for the
  null-model fit (used only when `method = "quasi"`).

## Value

A data frame with columns `logFC` (log2 fold change of `cluster` vs
`cluster_2`), `logCPM`, `Prop.1` and `Prop.2` (proportion of cells
expressing the gene in group 1 and group 2 respectively), `stat` (signed
quasi statistic), `PValue` (two-sided), and the adjusted p-value column.

## Examples

``` r

data("xenium_bc_spe")
spe <- normalizeAssay(spe)
spe <- runPCA(spe)
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
spe <- findNbrsSNN(spe, dimred = "PCA")
spe <- getClusters(spe, resolution = 0.5)
# cluster 1 vs cluster 2
de <- getDE(spe, cluster = 1, cluster_2 = 2)
# clusters 1 & 3 vs clusters 2 & 4
de2 <- getDE(spe, cluster = c(1, 3), cluster_2 = c(2, 4))
```
