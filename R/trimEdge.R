#' Trim spots/bins from the edges of a Visium (or VisiumHD) slide
#'
#' Removes whole lines of spots from the outer edges of an aligned Visium or
#' VisiumHD array, using the regular \code{array_row}/\code{array_col} grid.
#' Useful for quickly shaving off edge artefacts (folds, tears, capture
#' effects).
#'
#' @details
#' Edges are defined in the \strong{stored-coordinate} frame, which is
#' independent of any plotting choice: \code{left} = smallest x, \code{right} =
#' largest x, \code{top} = smallest y, \code{bottom} = largest y. This matches
#' the orientation of scider's default \emph{image} plot (which is Y-reversed).
#' Note that a plot drawn \emph{without} the Y-reversal (e.g. an object read
#' without an image, or \code{reverseY = FALSE}) is vertically mirrored, so
#' \code{top} then appears at the bottom of that plot.
#'
#' Trimming removes the \code{n} smallest/largest \code{array_col} values (for
#' \code{left}/\code{right}) and \code{array_row} values (for \code{top}/
#' \code{bottom}). Because of the Visium hexagonal offset, one \code{array_col}
#' value spans only alternate rows, so removing a single value trims roughly
#' half a visual column at a time; use \code{2} for a full column line.
#' \code{array_row} values are clean horizontal lines. VisiumHD's square grid
#' has no such half-offset.
#'
#' This function requires \code{array_row}/\code{array_col} in
#' \code{colData(spe)} and therefore only applies to Visium/VisiumHD.
#'
#' @param spe A SpatialExperiment with \code{array_row} and \code{array_col}
#'   in its \code{colData}.
#' @param trim Integer vector of length 4 giving the number of spot lines to
#'   remove from each edge, in the order \code{c(bottom, left, top, right)}
#'   (the base R \code{\link[graphics]{par}} \code{mar} order). Default
#'   \code{c(0, 0, 0, 0)} (no trimming).
#'
#' @return The SpatialExperiment with the edge spots removed.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' spe <- readVisium("path/to/visium/outs")
#' # 1 line off the bottom, 2 off the left, 3 off the top, 4 off the right
#' spe <- trimEdge(spe, trim = c(1, 2, 3, 4))
#' }
#'
trimEdge <- function(spe, trim = c(0, 0, 0, 0)) {
  if (is.null(spe$array_row) || is.null(spe$array_col)) {
    stop("trimEdge() requires 'array_row' and 'array_col' in colData(spe) ",
         "(Visium or VisiumHD).")
  }
  if (length(trim) != 4L) {
    stop("'trim' must be a length-4 vector: c(bottom, left, top, right).")
  }
  trim <- as.integer(trim)
  bottom <- trim[1]
  left   <- trim[2]
  top    <- trim[3]
  right  <- trim[4]

  ar <- spe$array_row
  ac <- spe$array_col
  coord <- spatialCoords(spe)
  x <- coord[, 1]
  y <- coord[, 2]

  # Map the array axes onto the stored-coordinate edges (left = min x,
  # right = max x, top = min y, bottom = max y), robust to how the indices
  # happen to be oriented relative to the coordinates.
  col_left_end <- if (stats::cor(ac, x) >= 0) "low" else "high"
  row_top_end  <- if (stats::cor(ar, y) >= 0) "low" else "high"
  opp <- function(end) if (end == "low") "high" else "low"

  # Return the n extreme (smallest or largest) distinct values of v to drop.
  pick <- function(v, n, end) {
    if (n <= 0L) return(v[0])
    u <- sort(unique(v))
    if (end == "low") head(u, n) else tail(u, n)
  }

  drop_col <- c(pick(ac, left,  col_left_end),
                pick(ac, right, opp(col_left_end)))
  drop_row <- c(pick(ar, top,    row_top_end),
                pick(ar, bottom, opp(row_top_end)))

  keep <- !(ac %in% drop_col | ar %in% drop_row)
  spe[, keep]
}
