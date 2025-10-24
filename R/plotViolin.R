#' Violin plot using genes or cell data
#' @param spe A SpatialExperiment object.
#' @param feature can be a vector of gene names in rownames(spe), or column 
#' names in colData(spe) if those columns are numeric.
#' @param assay Name of assay to use for plotting feature.
#' @param group.by values to group plot by. Must be in colData of spe and must
#' be either factor or character.
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
                       feature,
                       assay = "counts",
                       group.by = NULL,
                       type = c("raw","log","cpm","logcpm"),
                       point = FALSE,
                       color.by = NULL,
                       ncol = NULL,
                       pt.size = 0.3,
                       pt.alpha = 0.3,
                       pt.shape = ".",
                       label.y = "Expression") {
  # Retrieve feature
  dat <- lapply(feature, function(f) {
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
  dat <- data.frame(expression=unlist(dat), x=as.factor(rep(feature,each=ncol(spe))))
  
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
    group <- as.factor(rep(spe@colData[[group.by]],length(feature)))
  }
  
  ## Plotting
  p <- ggplot(data=dat) + geom_violin(aes(x=group,y=expression)) +
    labs(y=label.y, x="") + 
    theme_classic() +
    theme(axis.text.x=element_text(angle=-45,hjust=0)) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  # Separate by feature
  p <- p + facet_wrap(dat$x, ncol=ncol)
  
  # Separating by points/colors
  if (point) {
    # prepping color.by
    color <- NULL
    if (!is.null(color.by) && color.by %in% names(spe@colData)) {
      color <- as.factor(rep(spe@colData[[color.by]],length(feature)))
      }
    p <- p + geom_point(aes(x=group, y=expression, color=!!color),
                        position = position_jitter(seed=1, width=0.2),
                        shape = pt.shape,
                        size = pt.size,
                        alpha = pt.alpha) + 
      guides(color = guide_legend(override.aes = list(alpha=1,shape=19,size=1)))
  }
     
  p
}
