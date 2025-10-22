#' Violin plot using genes or cell data
#' @param spe A SpatialExperiment object.
#' @param y can be a gene name in rownames(spe) or cell data in coldata(spe)
#' @param assay Name of assay to use for plotting feature.
#' @param group.by values to group plot by. Must be in colData of spe.
#' @param type Transformation to apply for the group/feature. Options are "raw"
#' , "log", "cpm", "logcpm", or a function that accepts and returns a vector of 
#' the same length.
#' @param point Whether to plot points. 
#' @param color.by values to color points by. Must be in colData of spe.
#' @param ncol Number of column if group.by is used.
#' @param pt.size Size of points.
#' @param pt.alpha Alpha of points between 0 and 1.
#' @param pt.shape Shape of points.
#' @param label.y Label for the y-axis.
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' plotViolin(spe,c("cell_area","nucleus_area"),group.by="cell_type",label.y="Area")
plotViolin <- function(spe,
                       y,
                       assay = "counts",
                       group.by = NULL,
                       type = c("raw","log","cpm","logcpm"),
                       point = TRUE,
                       color.by = NULL,
                       ncol = NULL,
                       pt.size=0.3,
                       pt.alpha=0.3,
                       pt.shape=".",
                       label.y = "Expression") {
  # Retrieve y
  dat <- lapply(y, function(f) {
    if (f %in% rownames(spe)) {
      SummarizedExperiment::assay(spe,assay)[f,]
    } else if (f %in% names(spe@colData)) {
      spe@colData[[f]]
    } else {
      message(paste0("Couldn't find ",f,". Skipping."))
      NULL
    }
  })
  
  # transform dat to long matrix
  dat <- data.frame(expression=unlist(dat),x=as.factor(rep(y,each=ncol(spe))))
  
  # transform expression
  if (is.character(type)) {
    type <- switch(match.arg(type),
                   raw = NULL,
                   log = function(x) {log2(x+1)},
                   cpm = function(x) {
                     (x+0.5)/colSums(as.matrix(spe@assays@data[[assay]]))*1e6
                   },
                   logcpm = function(x) {
                     log2((x+0.5)/colSums(as.matrix(spe@assays@data[[assay]]))*1e6)
                   })
  }
  if(is.function(type)) {
    tryCatch({dat$expression <- type(dat$expression)},
             error = function(e){
               message("Error when applying 'type'. Skipping 'type'.")
             })
  }

  group <- ""
  if (!is.null(group.by) && group.by %in% names(spe@colData)) {
    group <- as.factor(rep(spe@colData[[group.by]],length(y)))
  }
  
  ## Plotting
  p <- ggplot(data=dat) + geom_violin(aes(x=group,y=expression)) +
    theme(axis.text.x=element_text(angle=-45,hjust=0)) + 
    labs(y=label.y,
         x="")
  
  # Separate by feature
  p <- p + facet_wrap(dat$x,ncol=ncol)
  
  # Separating by points/colors
  if (point) {
    # prepping color.by
    color <- NULL
    if (!is.null(color.by) && color.by %in% names(spe@colData)) {
      color <- as.factor(rep(spe@colData[[color.by]],length(y)))
      }
    p <- p + geom_point(aes(x=group,y=expression,color=!!color),
                        position = position_jitter(seed = 1,width=0.2),
                        shape = pt.shape,
                        size = pt.size,
                        alpha = pt.alpha) + 
      guides(color = guide_legend(override.aes = list(alpha=1,shape=19,size=1)))
  }
     
  p
}
