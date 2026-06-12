rownames2col <- function(df, rn = "rowname") {
    stopifnot(is.data.frame(df))
    df[[rn]] <- rownames(df)
    rownames(df) <- NULL
    return(df)
}

col2rownames <- function(df, rn = "rowname") {
    stopifnot(is.data.frame(df))
    df <- as.data.frame(df)
    rownames(df) <- df[[rn]]
    df[[rn]] <- NULL
    return(df)
}

col.spec <- c(
    "#D53E4F", "#F46D43", "#FDAE61",
    "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4",
    "#66C2A5", "#3288BD", "#5E4FA2"
)

col.pMedium <- c(
    "#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9",
    "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA"
)

col.pDark <- c(
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
)

col.pLight <- c(
    "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
    "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5"
)

col.p10 <- col.pMedium
col.p20 <- c(col.pDark, col.pLight)
col.p30 <- c(col.pDark, col.pLight, col.pMedium)

selectColor <- function(n) {
    if (n <= 10) {
        return(col.p10[seq_len(n)])
    } else if (n <= 20) {
        return(col.p20[seq_len(n)])
    } else if (n <= 30) {
        return(col.p30[seq_len(n)])
    } else {
        return(rep_len(col.p30, n))
    }
}

# Build a colour palette for plotting points by group.
# For discrete groups, the special level "unassigned" is coloured black and
# excluded from the palette, so the real clusters keep contiguous palette
# colours. Returns an unnamed vector for continuous groups and a named vector
# (keyed by level) for discrete groups.
.buildColP <- function(group, cols, isContinuous,
                       special = "unassigned", special_col = "black") {
  if (isContinuous) {
    n <- length(unique(group))
    if (is.null(cols)) return(col.spec)
    if (is.function(cols)) return(cols(n))
    return(cols)
  }
  lv <- levels(as.factor(group))
  has_special <- special %in% lv
  main_lv <- if (has_special) setdiff(lv, special) else lv

  if (!is.null(cols) && !is.function(cols) && !is.null(names(cols))) {
    # Named colour vector: map colours to levels by name. Levels not named fall
    # back to the default palette (or black for 'unassigned').
    col.p <- rep(NA_character_, length(lv))
    names(col.p) <- lv
    matched <- intersect(lv, names(cols))
    col.p[matched] <- cols[matched]
    unmatched <- main_lv[!main_lv %in% matched]
    if (length(unmatched)) col.p[unmatched] <- selectColor(length(unmatched))
  } else {
    # Unnamed vector, function, or NULL: assign palette to the retained levels
    # in level order.
    col.p <- if (is.null(cols)) selectColor(length(main_lv))
             else if (is.function(cols)) cols(length(main_lv))
             else rep_len(cols, length(main_lv))
    names(col.p) <- main_lv
  }
  # 'unassigned' is black by default unless the user supplied a colour for it.
  if (has_special && is.na(col.p[special])) col.p[special] <- special_col
  col.p
}

# Prepare colour, size and draw-order for highlighting a subset of cells.
# `highlight` is either (a) a vector of group.by levels (characters or cluster
# numbers) to emphasise, or (b) a logical vector of length ncells selecting cells
# directly (e.g. counts(spe)["Sox9", ] >= 3). Non-highlighted cells are drawn
# first in light grey at `pt.size`; highlighted cells are drawn last (shuffled
# among themselves, so no level sits systematically on top) at
# `pt.size.highlight`. Highlight colours come from `cols.highlight` (a single
# colour for all, a vector of one colour per 'highlight' entry, or NULL to keep
# the usual group.by colours). Returns the re-levelled grouping factor, the
# matching palette, a per-cell size vector (all in the original cell order) and
# the row draw order.
.prepHighlight <- function(group, highlight, cols, cols.highlight = NULL,
                           pt.size, pt.size.highlight, grey = "grey80") {
  if (!is.null(group)) group <- as.factor(group)

  if (is.logical(highlight)) {
    # Logical mask: select cells directly. A single highlight category.
    if (!is.null(group) && length(highlight) != length(group))
      stop("Logical 'highlight' must have length equal to the number of cells (",
           length(group), ").")
    is_hl   <- highlight & !is.na(highlight)
    hcat    <- "highlight"
    hl_cats <- "highlight"
    hl_req  <- "highlight"
    # A mask has no group level, so there is no "usual" palette colour for it.
    if (is.null(cols.highlight)) cols.highlight <- "red"
  } else {
    if (is.null(group))
      stop("'highlight' given as levels requires 'group.by' to be set.")
    hl_req  <- unique(as.character(highlight))   # user-supplied order
    is_hl   <- as.character(group) %in% hl_req
    hcat    <- as.character(group)
    hl_cats <- hl_req[hl_req %in% levels(group)]  # valid levels, user order
  }
  if (!any(is_hl)) warning("'highlight' matched no cells.")

  # Highlight colours, from cols.highlight:
  #  - NULL          -> the usual group.by palette colours
  #  - single colour -> all highlighted cells that colour (default "red")
  #  - vector        -> one colour per 'highlight' entry (matched by position)
  if (is.null(cols.highlight)) {
    base_p <- .buildColP(group, cols, isContinuous = FALSE)
    hl_col <- base_p[hl_cats]
  } else if (length(cols.highlight) == 1) {
    hl_col <- rep(cols.highlight, length(hl_cats))
  } else {
    if (length(cols.highlight) != length(hl_req))
      stop("'cols.highlight' must be length 1 or match the number of 'highlight' ",
           "entries (", length(hl_req), ").")
    cmap <- cols.highlight
    names(cmap) <- hl_req
    hl_col <- cmap[hl_cats]
  }

  other_lab <- "other"
  while (other_lab %in% c(hl_cats, levels(group))) other_lab <- paste0(other_lab, ".")

  new_chr <- ifelse(is_hl, hcat, other_lab)
  group2  <- factor(new_chr, levels = c(hl_cats, other_lab))

  col.p <- character(0)
  col.p[hl_cats]   <- hl_col
  col.p[other_lab] <- grey

  size <- ifelse(is_hl, pt.size.highlight, pt.size)

  # Draw non-highlighted first, then highlighted last (shuffled among themselves).
  hl_idx <- which(is_hl)
  if (length(hl_idx) > 1) hl_idx <- sample(hl_idx)
  ord <- c(which(!is_hl), hl_idx)

  list(group = group2, col.p = col.p, size = size, order = ord)
}

# Order cluster factor levels: numeric labels first (numerically), then any
# non-numeric labels (e.g. merged "2&5&7") alphabetically, with `special`
# ("unassigned") always last.
.orderClusterLevels <- function(labels, special = "unassigned") {
  u <- setdiff(unique(labels), special)
  num <- suppressWarnings(as.numeric(u))
  is_num <- !is.na(num)
  c(u[is_num][order(num[is_num])],
    sort(u[!is_num]),
    if (special %in% labels) special)
}

col.lisa <- c("#eeeeee", "#FF0000", "#0000FF", "#a7adf9",
              "#f4ada8", "#464646", "#999999")
col.pval <- c("#3644E5", "#FFFFBF", "#FF5D53")

contour_brks <- getFromNamespace("contour_breaks", "ggplot2")
unique00 <- getFromNamespace("unique0", "ggplot2")
data_frame00 <- getFromNamespace("data_frame0", "ggplot2")

colsum <- function(x, group, reorder = TRUE,...) {
  return(Matrix::t(rowsum(Matrix::t(x), group, reorder = reorder, ...)))
}

# Clean up vector of names. If "overall" is in names (after cleaning up) then 
# also removes "overall" from it if more than 1 names.
cleanName <- function(names) {
  names <- janitor::make_clean_names(names)
  names <- gsub("^density_|_contour$|_roi$","",names)
  if ("overall" %in% names && length(names)>1) names <- names[names!="overall"]
  names <- sort(names)
  return(names)
}