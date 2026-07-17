#' Dot plot of gene expression by groups
#' 
#' Visualizing average expression and percentage of cell that expression a 
#' certain gene(s). Similar to Seurat's DotPlot.
#' @param spe A SpatialExperiment object.
#' @param feature vector of feature names.
#' @param assay Name of assay to use for plotting feature.
#' @param group.by values to group plot by. Must be in colData of spe.
#' @param detection.limit threshold for minimum expression value for percentage 
#' expression calculation (dot size)
#' @param range Upper and lower bound for average expression. Values
#' beyond this range are snapped to this limit. If only one value is provided,
#' it is taken as the upper bound.
#' @param scale Whether to scale the average expression data of each feature
#' using \link[base]{scale}.
#' @param cols Custom colour palette (a vector of colours defining the fill
#' gradient). Default NULL uses a light-grey (low) to dark-red (high) gradient.
#' @param dot.scale scale the radius of the plot. See \link[ggplot2]{scale_radius}
#' @param flip.axes Whether to flip the axes.
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' plotDots(spe,feature = rownames(spe)[1:5])
plotDots <- function(spe,
                     feature,
                     assay = "counts",
                     group.by = "cell_type",
                     detection.limit = 0, # Dot size
                     range = c(-Inf,Inf), # Dot colour
                     scale = TRUE,
                     cols = NULL,
                     dot.scale = 6,
                     flip.axes = FALSE) {
  # Get features
  f_missing <- !(feature %in% rownames(spe))
  if (any(f_missing)) {
    message(paste0(paste(feature[f_missing],collapse=", "),
                   " not found. Skipping"))
    feature <- feature[!f_missing]
  }
  exprs <- SummarizedExperiment::assay(spe,assay)[feature,,drop=FALSE]
  # Get group
  if (!is.null(group.by) && group.by %in% names(spe@colData)) {
    group <- spe[[group.by]]
  } else {
    group <- ""
  }
  group <- as.factor(group)
  
  # Calculate mean & percentage expression of each group.
  dat <- stats::aggregate.data.frame(as.matrix(Matrix::t(exprs)), list(group),
                                     FUN=function(x){
                                       n=length(x)
                                       mean = log1p(mean(expm1(x)))
                                       percentage = sum(x>detection.limit)/n
                                       matrix(c(mean, percentage),
                                              ncol=2)})
  if (scale) {
    for (i in 2:(length(feature)+1)) {
      dat[[i]][,1] = as.vector(scale(dat[[i]][,1]))
    }
  }
  
  # Convert to long
  dat <- data.frame(rep(dat[[1]], ncol(dat)-1), # Group
                    rep(colnames(dat)[-1], each = nrow(dat)), # Feature
                    do.call(rbind,dat[,-1,drop=FALSE])) # Summary stats
  colnames(dat) = c(group.by,"feature","average","percentage")
  
  # Threshold average expression
  if (length(range)==1) {range = c(-Inf,range)}
  dat$average[dat$average<range[1]] <- range[1]
  dat$average[dat$average>range[2]] <- range[2]
  
  # Convert to factor to keep the specified order when plotted
  dat$feature <- factor(dat$feature,levels = feature)
  
  # Plotting
  p <- ggplot(data=dat) +
    geom_point(aes(.data[[group.by]],feature,size=.data[["percentage"]],fill=.data[["average"]]),
               shape=21, colour="black", stroke=0.3) +
    theme_minimal() +
    theme(axis.line=element_line(colour="black"),
          axis.text.x=element_text(angle=45,hjust=1)) +
    scale_radius(range=c(0,dot.scale))
  if (is.null(cols)) {
    p <- p + scale_fill_gradient(low="lightgrey", high="darkred")
  } else {
    p <- p + scale_fill_gradientn(colours=cols)
  }
  if (flip.axes) {
    p <- p + coord_flip()
  }
  
  return(p)
}