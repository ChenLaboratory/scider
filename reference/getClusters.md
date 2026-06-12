# Cluster cells in spe using graph methods.

Cluster cells in spe using graph methods.

## Usage

``` r
getClusters(
  spe,
  nbrs_name = NULL,
  method = c("leiden", "louvain"),
  resolution = 1,
  cluster_name = "cluster",
  seed = 1,
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- nbrs_name:

  Name of neighbour list for clustering. If NULL, will use the newest
  one in spe@metadata\$nbrs\$cell or create one if none are available.

- method:

  Clustering methods. Options are leiden and louvain.

- resolution:

  Higher resolution for more clusters and lower for fewer clusters. See
  [cluster_leiden](https://r.igraph.org/reference/cluster_leiden.html)
  and
  [cluster_louvain](https://r.igraph.org/reference/cluster_louvain.html)

- cluster_name:

  Name to store the clusters in spe's
  [colData](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

- seed:

  seed for clustering

- ...:

  Other clustering arguments for
  [cluster_leiden](https://r.igraph.org/reference/cluster_leiden.html)
  or
  [cluster_louvain](https://r.igraph.org/reference/cluster_louvain.html)

## Value

A spe with the clusters stored in
[reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html).

A SpatialExperiment object

## Details

Cluster cells with igraph using SNN calculated by
[findNbrsSNN](https://chenlaboratory.github.io/scider/reference/findNbrsSNN.md).
Any neighbour list in spe@metadata\$nbrs\$cell can also be used

## Examples

``` r

data("xenium_bc_spe")
spe <- normalizeAssay(spe)
spe <- runPCA(spe)
#> Genes with 0 variance are excluded: ENSG00000135218 NegControlProbe_00002 NegControlCodeword_0504 NegControlCodeword_0509 NegControlCodeword_0510 NegControlCodeword_0511 NegControlCodeword_0512 NegControlCodeword_0516 NegControlCodeword_0517 NegControlCodeword_0518 NegControlCodeword_0519 NegControlCodeword_0520 NegControlCodeword_0522 NegControlCodeword_0526 NegControlCodeword_0527 NegControlCodeword_0530 NegControlCodeword_0536 NegControlCodeword_0537 BLANK_0030 BLANK_0163 BLANK_0165 BLANK_0212 BLANK_0221 BLANK_0230 BLANK_0237 BLANK_0311 BLANK_0361 BLANK_0365 BLANK_0382 BLANK_0384 BLANK_0387 BLANK_0388 BLANK_0391 BLANK_0393 BLANK_0397 BLANK_0399 BLANK_0404 BLANK_0406 BLANK_0410 BLANK_0411 BLANK_0418 BLANK_0425 BLANK_0432 BLANK_0447
spe <- findNbrsSNN(spe,dimred="PCA")
#> [1] "Getting K-nearest neighbour"
#> [1] "Getting shared nearest neighbour"
spe <- getClusters(spe, resolution=0.5)
#> 1 cells in clusters with <= 1 cells labelled 'unassigned'.
```
