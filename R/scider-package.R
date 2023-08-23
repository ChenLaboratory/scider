#' @import ggplot2
#' @importFrom methods is
#' @import shiny
#' @importFrom plotly plot_ly
#' @importFrom plotly renderPlotly
#' @importFrom plotly add_markers
#' @importFrom plotly layout
#' @importFrom plotly event_data
#' @importFrom SummarizedExperiment colData
#' @import SpatialExperiment

NULL


#' Spatial cell-type inter-correlation by density in R.
#'
#' `scider` implements functions to analyse spatial transcriptomics data with 
#' cell type annotations by performing cell type correlation via density estimation
#' and cell type co-localization via real number distance. Functions include 
#' density estimation, statistical modelling and visualizations.
#'
#' `scider` uses SpatialExperiment objects as the main infrastructure, which can
#' easily be integrated with a wide variety of Bioconductor packages.
#'
#'
#' @author Ning Liu \email{liu.n@@wehi.edu.au}, Mengbo Li \email{li.me@@wehi.edu.au},
#' Yunshun Chen \email{yuchen@@wehi.edu.au}
#' @name scider-package
#' @docType package
#' @aliases scider scider-package
#' @keywords internal
#'
NULL


#' Description of the scider example datasets
#'
#' scider-package has 1 datasets: \itemize{
#'   \item xenium_bc_spe Example subset of a Xenium breast cancer dataset,
#'  }
#'
#' @docType data
#' @name xenium_bc_spe
#' @usage data("xenium_bc_spe")
#' @keywords internal
#' @format A SpatialExperiment object
#' @examples
#' data(xenium_bc_spe)
"spe"