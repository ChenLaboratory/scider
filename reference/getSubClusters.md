# Sub-cluster a single cluster.

Take one cluster from an existing clustering and partition it further by
community detection on the sub-graph of its cells, similar to
`Seurat::FindSubCluster`. New sub-clusters are labelled `<cluster>_1`,
`<cluster>_2`, ... (largest first).

## Usage

``` r
getSubClusters(
  spe,
  cluster,
  resolution = 1,
  nbrs_name = NULL,
  method = c("leiden", "louvain"),
  cluster_name = "cluster",
  new_name = cluster_name,
  relabel = FALSE,
  min_size = NULL,
  start_from = NULL,
  sep = "_",
  seed = 1,
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- cluster:

  The cluster label to sub-cluster (e.g. 2 or "2").

- resolution:

  Resolution for the sub-clustering. Higher gives more sub-clusters. See
  [cluster_leiden](https://r.igraph.org/reference/cluster_leiden.html)
  and
  [cluster_louvain](https://r.igraph.org/reference/cluster_louvain.html).

- nbrs_name:

  Name of neighbour list to use. If NULL, uses the newest one in
  spe@metadata\$nbrs\$cell.
  [findNbrsSNN](https://chenlaboratory.github.io/scider/reference/findNbrsSNN.md)
  must have been run first.

- method:

  Clustering method. Options are leiden and louvain.

- cluster_name:

  Name of the cluster column in
  [colData](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  to sub-cluster.

- new_name:

  Name of the column to store the result. Defaults to overwriting
  cluster_name.

- relabel:

  Logical. If TRUE, renumber all clusters and sub-clusters by size
  (largest first) after sub-clustering, replacing the `<cluster>_<n>`
  labels with plain numbers. "unassigned" is kept as the last level.
  Defaults to FALSE.

- min_size:

  Sub-clusters with at most this many cells are merged into the larger
  sub-clusters (graph connection first, then nearest centroid), as in
  [getClusters](https://chenlaboratory.github.io/scider/reference/getClusters.md).
  Defaults to NULL, which uses either 5 or 0.01\\ the cells, whichever
  is smaller.

- start_from:

  Integer at which numbering starts, for both the sub-cluster suffixes
  and the relabel renumbering. Defaults to NULL, which infers the base
  (0- or 1-based) from the existing labels.

- sep:

  Separator between a cluster and its sub-cluster number. Default "\_".

- seed:

  Seed for clustering.

- ...:

  Other clustering arguments for
  [cluster_leiden](https://r.igraph.org/reference/cluster_leiden.html)
  or
  [cluster_louvain](https://r.igraph.org/reference/cluster_louvain.html).

## Value

A SpatialExperiment with the sub-clusters stored in colData.

## Details

The sub-graph of the target cluster's cells is induced from the SNN
graph and clustered on its own. Cells that are weakly connected within
the cluster may form small or singleton sub-clusters; lower the
resolution if this is not wanted.

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
#> Reassigned all 1 cells from clusters with <= 1 cells (0 by graph, 1 by distance).
# Split cluster 2 into 2_1, 2_2, ...
spe <- getSubClusters(spe, cluster = 2, resolution = 0.5)
#> Split cluster '2' into 4 sub-clusters.
```
