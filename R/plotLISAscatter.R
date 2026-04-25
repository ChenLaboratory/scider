#' Scatterplot for local moran's I
#'
#' plot result obtained from localMoran()
#' @param lisa A list obtained from \link[sciderHex]{localMoran}
#' @param quadrant.count Whether to count values at each quadrant (Low-Low, 
#' Low-High, High-High, High-Low) 
#' @param text.size Numeric for text size of quadrant.count
#' @param xlab label for the x-axis
#' @param ylab label for the y-axis
#' 
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' data("xenium_bc_spe")
#' dat <- spe$total_counts
#' spe <- findNbrsSpatial(spe,k=50)
#' res <- localMoran(spe,data1=dat,at="cell")
#' plotLISAscatter(res)
plotLISAscatter <- function(lisa, 
                            quadrant.count = TRUE,
                            text.size=11/.pt,
                            xlab="Data",
                            ylab="Spatial Lag") {
  if (is.null(lisa$spatiallag)) {
    stop("Missing spatiallag. LISA appears to be global instead of local")
  }
  
  keep <- lisa$cluster!="Undefined"
  dat <- (lisa$lisa/lisa$spatiallag)[keep]
  lag <- lisa$spatiallag[keep]
  Cluster <- lisa$cluster[keep]
  cols <- col.lisa[sort(unique(Cluster))]
  globalMoran <- mean(lisa$lisa)
  
  p <- ggplot2::ggplot() + 
    geom_point(aes(x=dat,y=lag,color = Cluster)) + 
    scale_color_manual(values = cols) +
    geom_abline(intercept=0,slope=globalMoran) +
    geom_hline(yintercept=0,linetype="dotted") +
    geom_vline(xintercept=0,linetype="dotted") +
    ggtitle(sprintf("Moran's I: %0.3f",globalMoran))+
      # paste("Moran's I:",globalMoran)) +
    labs(x = xlab, y = ylab) +
    theme_classic()
  
  # Count for each quadrants
  if (quadrant.count) {
    p <- p +
      # Going clockwise from bottom left 
      geom_text(
        size = text.size,
        aes(
          x=c(-Inf,-Inf,Inf,Inf),
          y=c(-Inf,Inf,Inf,-Inf),
          label=c(paste0("Low-Low (",sum(dat<0&lag<0),")"),
                  paste0("Low-High (",sum(dat<0&lag>0),")"),
                  paste0("High-High (",sum(dat>0&lag>0),")"),
                  paste0("High-Low (",sum(dat>0&lag<0),")")
          ),
          hjust=c(0,0,1,1),
          vjust=c(-0.2,1,1,-0.2)
        )
      )
  }
  p
}
