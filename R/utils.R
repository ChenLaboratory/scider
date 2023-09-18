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
