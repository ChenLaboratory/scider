# Dot plot or heatmap of the top marker genes per cluster.

Takes the output of
[getMarkers](https://chenlaboratory.github.io/scider/reference/getMarkers.md)
and visualises the top `n` marker genes of each cluster, either as a dot
plot or a genes-by-clusters expression heatmap. Genes are block-ordered
by cluster (in the order of `names(markers)`); a gene that is a top
marker for more than one cluster is shown once, under the first cluster
it appears in.

## Usage

``` r
plotTopMarkers(
  markers,
  spe,
  n = 5,
  min.prop = 0,
  type = c("dotplot", "heatmap"),
  profile = c("pseudobulk", "cell"),
  cluster_name = "cluster",
  assay = "counts",
  scale = TRUE,
  range = NULL,
  cols = NULL,
  dot.scale = 6,
  ...
)
```

## Arguments

- markers:

  The named list returned by
  [getMarkers](https://chenlaboratory.github.io/scider/reference/getMarkers.md)
  (one data frame per cluster, gene names as row names).

- spe:

  The SpatialExperiment object the markers were computed from.

- n:

  Integer. Number of top markers to take from each cluster. Default 5.

- min.prop:

  Numeric in \[0, 1\]. Drop candidate markers whose detection rate
  (proportion of cells expressing) in their cluster is below this value
  before taking the top `n`, so low-prevalence genes are skipped and
  back-filled from further down the ranking. Default 0 (no filtering).

- type:

  Character. "dotplot" (default) or "heatmap".

- profile:

  Character. "pseudobulk" (default) or "cell"; see Details.

- cluster_name:

  Character. Column in `colData(spe)` holding the cluster labels used in
  `markers`. Default "cluster".

- assay:

  Name of the assay to use for expression. Only used when
  `profile = "cell"` (pseudo-bulk always uses raw counts). Default
  "counts".

- scale:

  Logical. Scale expression per gene (z-score across clusters for the
  colour). Default TRUE.

- range:

  Numeric. Lower and upper caps for the colour values; values outside
  are clamped (winsorised). A single value is taken as the upper bound.
  Default NULL picks `c(-1.5, 2.2)` when `scale = TRUE` (suited to
  z-scores) and `c(-Inf, Inf)` (no capping) otherwise.

- cols:

  Custom colour palette. For a cell-level dot plot a gradient passed to
  [plotDots](https://chenlaboratory.github.io/scider/reference/plotDots.md);
  otherwise the full colour vector for the gradient (dot plot) or
  [pheatmap](https://rdrr.io/pkg/pheatmap/man/pheatmap.html) (heatmap).

- dot.scale:

  Numeric. Scales the dot radius. See
  [scale_radius](https://ggplot2.tidyverse.org/reference/scale_size.html).
  Default 6.

- ...:

  Further arguments passed to
  [plotDots](https://chenlaboratory.github.io/scider/reference/plotDots.md)
  (cell-level dot plot) or
  [pheatmap](https://rdrr.io/pkg/pheatmap/man/pheatmap.html) (heatmap).

## Value

A ggplot object when `type = "dotplot"`, or a pheatmap object when
`type = "heatmap"`.

## Details

Two expression profiles are available via `profile`:

- **"pseudobulk"** (default): per-cluster log-CPM from
  [spe2PB](https://chenlaboratory.github.io/scider/reference/spe2PB.md) +
  [cpm](https://rdrr.io/pkg/edgeR/man/cpm.html) - the same space the
  markers were tested in. Dot colour is (scaled) log-CPM; dot size is
  the proportion of cells expressing the gene.

- **"cell"**: Seurat-style cell-level summaries via
  [plotDots](https://chenlaboratory.github.io/scider/reference/plotDots.md)
  (dot plot) or cell-level mean expression (heatmap).

## Examples

``` r

data("xenium_bc_spe")
spe <- normalizeAssay(spe)
spe <- runPCA(spe)
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
spe <- findNbrsSNN(spe, dimred = "PCA")
spe <- getClusters(spe, resolution = 0.5)
markers <- getMarkers(spe)
plotTopMarkers(markers, spe, n = 5)

plotTopMarkers(markers, spe, n = 5, type = "heatmap")
```
