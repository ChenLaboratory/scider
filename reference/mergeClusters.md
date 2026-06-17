# Manually merge clusters.

Combine specified clusters into a single cluster, optionally relabelling
all clusters afterwards. This is a manual alternative to the
modularity-based mergeCommunities: the user chooses exactly which
clusters to combine.

## Usage

``` r
mergeClusters(
  spe,
  merge,
  cluster_name = "cluster",
  new_name = cluster_name,
  merge_label = NULL,
  relabel = FALSE,
  start_from = NULL,
  sep = "&"
)
```

## Arguments

- spe:

  A SpatialExperiment object.

- merge:

  A vector of cluster labels to merge into one (e.g. c(2, 5, 7)), or a
  list of such vectors to perform several merges in one call.

- cluster_name:

  Name of the cluster column in
  [colData](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  to merge.

- new_name:

  Name of the column to store the result. Defaults to overwriting
  cluster_name.

- merge_label:

  Optional label for the merged cluster (single merge only). Defaults to
  the merged cluster labels joined by `sep`, e.g. "2&5&7".

- relabel:

  Logical. If TRUE, renumber all clusters by size (largest first) after
  merging, replacing the merged label with a number. "unassigned" is
  always kept as the last level. Defaults to FALSE (keep the merged
  label).

- start_from:

  Integer at which renumbering starts when relabel=TRUE. Defaults to
  NULL, which infers the base (0- or 1-based) from the existing labels,
  matching what getClusters produced. Ignored when relabel=FALSE.

- sep:

  Separator for the auto-generated merged label. Default "&".

## Value

A SpatialExperiment with the merged clusters stored in colData.

## Details

Clusters listed in `merge` that are not present in cluster_name are
ignored; a merge group with fewer than two present clusters is skipped
with a warning. Any "unassigned" cells are left untouched (and kept as
the last factor level) unless explicitly included in `merge`.

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
# Merge clusters 2, 5 and 7 into one labelled "2&5&7"
spe <- mergeClusters(spe, merge = c(2, 5, 7))
# Merge and then renumber everything 1..K by size
spe <- mergeClusters(spe, merge = c(2, 5, 7), relabel = TRUE)
#> Warning: Merge group {2, 5, 7} has fewer than two clusters present; skipping.
```
