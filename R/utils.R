rownames2col <- function(df, rn = "rowname"){
  stopifnot(is.data.frame(df))
  df[[rn]] <- rownames(df)
  rownames(df) <- NULL
  return(df)
}

col2rownames <- function(df, rn = "rowname"){
  stopifnot(is.data.frame(df))
  df <- as.data.frame(df)
  rownames(df) <- df[[rn]]
  df[[rn]] <- NULL
  return(df)
}


col.spec <- c("#D53E4F", "#F46D43", "#FDAE61",
              "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4",
              "#66C2A5", "#3288BD", "#5E4FA2")

