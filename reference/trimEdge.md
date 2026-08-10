# Trim spots/bins from the edges of a Visium (or VisiumHD) slide

Removes whole lines of spots from the outer edges of an aligned Visium
or VisiumHD array, using the regular `array_row`/`array_col` grid.
Useful for quickly shaving off edge artefacts (folds, tears, capture
effects).

## Usage

``` r
trimEdge(spe, trim = c(0, 0, 0, 0))
```

## Arguments

- spe:

  A SpatialExperiment with `array_row` and `array_col` in its `colData`.

- trim:

  Integer vector of length 4 giving the number of spot lines to remove
  from each edge, in the order `c(bottom, left, top, right)` (the base R
  [`par`](https://rdrr.io/r/graphics/par.html) `mar` order). Default
  `c(0, 0, 0, 0)` (no trimming).

## Value

The SpatialExperiment with the edge spots removed.

## Details

Edges are defined in the **stored-coordinate** frame, which is
independent of any plotting choice: `left` = smallest x, `right` =
largest x, `top` = smallest y, `bottom` = largest y. This matches the
orientation of scider's default *image* plot (which is Y-reversed). Note
that a plot drawn *without* the Y-reversal (e.g. an object read without
an image, or `reverseY = FALSE`) is vertically mirrored, so `top` then
appears at the bottom of that plot.

Trimming removes the `n` smallest/largest `array_col` values (for
`left`/`right`) and `array_row` values (for `top`/ `bottom`). Because of
the Visium hexagonal offset, one `array_col` value spans only alternate
rows, so removing a single value trims roughly half a visual column at a
time; use `2` for a full column line. `array_row` values are clean
horizontal lines. VisiumHD's square grid has no such half-offset.

This function requires `array_row`/`array_col` in `colData(spe)` and
therefore only applies to Visium/VisiumHD.

## Examples

``` r

if (FALSE) { # \dontrun{
spe <- readVisium("path/to/visium/outs")
# 1 line off the bottom, 2 off the left, 3 off the top, 4 off the right
spe <- trimEdge(spe, trim = c(1, 2, 3, 4))
} # }
```
