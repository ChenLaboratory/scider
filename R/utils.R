rownames2col <- function(df, rn = "rwoname"){
  stopifnot(is.data.frame(df))
  df[[rn]] <- rownames(df)
  rownames(df) <- NULL
  return(df)
}
