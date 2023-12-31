% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotContourRegion.R
\name{plotContourRegion}
\alias{plotContourRegion}
\title{Visualising an sf object (for internal use only at the moment)}
\usage{
plotContourRegion(
  spe,
  coi,
  id = "cell_type",
  overlay = c("density", "cell"),
  sub.level
)
}
\arguments{
\item{spe}{A SpatialExperiment object.}

\item{coi}{A character vector of length 1 of the cell type of interest.}

\item{id}{A character. The name of the column of colData(spe) containing
the cell type identifiers. Set to cell_type by default.}

\item{overlay}{Character vector. Either plot overlay on density or cells.}

\item{sub.level}{Numeric vector of length 1 or 2, identifies which
density level to plot. When length is 1, plot the density region above
this level. When length is 2, plot the density region
between the two levels.}
}
\value{
A ggplot object.
}
\description{
Visualising an sf object (for internal use only at the moment)
}
