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
#' @param cols Colour palette for violins. Can be a vector of colours or a 
#' function that accepts an integer n and return n colours.
#' @param ncol Number of column if group.by is used.
#' @param pt.size Size of points.
#' @param pt.alpha Alpha of points between 0 and 1.
#' @param pt.shape Shape of points.
#' @param ylab Label for the y-axis.
#' @param xlab Label for the x-axis
#' @export
#' @examples
#'
#' data("xenium_bc_spe")
#' plotViolin(spe,c("cell_area","nucleus_area"),group.by="cell_type",ylab="Area")
plotViolin <- function(spe,
                       feature,
                       assay = "counts",
                       group.by = NULL,
                       type = c("raw","log","cpm","logcpm"),
                       point = FALSE,
                       cols = NULL,
                       ncol = NULL,
                       pt.size = 0.3,
                       pt.alpha = 0.3,
                       pt.shape = ".",
                       ylab = "Expression",
                       xlab = NULL) {
  # Retrieve feature
  dat = list()
  for (f in feature) {
    if (f %in% rownames(spe)) {
      d <- SummarizedExperiment::assay(spe,assay)[f,]
    } else if (f %in% names(spe@colData)) {
      d <- spe@colData[[f]]
    } else {
      message(paste0("Couldn't find ",f,". Skipping."))
      d <- NULL
    }
    # Check for non-numeric
    if (!is.null(d) && !is.numeric(d)) {
      message(paste0(f," is non-numeric. Skipping."))
      d <- NULL
    }
    dat[[f]] <- d
  }
  if (length(dat) == 0) stop("No valid feature found")
  feature <- names(dat)
  
  # Transform dat to long matrix
  dat <- data.frame(expression=unlist(dat), x=factor(rep(feature,each=ncol(spe)),levels = feature))
  
  # Group to separate feature by
  if (!is.null(group.by)) {
    if (!group.by %in% names(spe@colData)) {
      stop(sprintf("Couldn't find %s in colData",group.by))
    } 
    if (is.numeric(spe@colData[[group.by]])) {
      stop(paste0(group.by," must be either factor or character."))
    }
    dat$group <- as.factor(rep(spe@colData[[group.by]],length(feature)))
  } else {
    dat$group <- ""
  }
  
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
  
  # Colours
  n_color <- length(unique(dat$group))
  if (is.null(cols)) col.p <- selectColor(n_color)
  else if (is.function(cols)) {
    col.p <- cols(n_color)
  } else {# cols is vector
    col.p <- rep_len(cols,n_color)
  }

  ## Plotting
  p <- ggplot(data=dat, aes(x=.data[["group"]],y=.data[["expression"]])) + 
    geom_violin(aes(fill=.data[["group"]])) +
    labs(y=ylab, x=xlab) + 
    theme_classic() +
    theme(axis.text.x=element_text(angle=-45,hjust=0),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none") + 
    scale_fill_manual(values=col.p)

  # Separate by feature
  p <- p + facet_wrap(dat$x, ncol=ncol)
  
  # Separating by points/colours
  if (point) {
    p <- p + geom_point(position = position_jitter(seed=1, width=0.2),
                        shape = pt.shape,
                        size = pt.size,
                        alpha = pt.alpha)
  }
  p
}
