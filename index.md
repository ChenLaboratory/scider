# scider: Spatial cell-type inter-correlation by density in R. ![](reference/figures/scider_sticker.png)

[![R-CMD-check](https://github.com/ChenLaboratory/scider/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/ChenLaboratory/scider/actions)
[![Codecov test
coverage](https://codecov.io/gh/ChenLaboratory/scider/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/ChenLaboratory/scider?branch=devel)

## Overview

*scider* implements functions to analyse spatial transcriptomics data
with cell type annotations by performing cell type correlation via
density estimation and cell type co-localization via real number
distance. Functions include density estimation, statistical modelling
and visualizations.

## Installation

Install the released version from Bioconductor:

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scider")
```

Install the development version from GitHub:

``` r

library(devtools)
devtools::install_github("ChenLaboratory/scider")
```

## Case studies

Step-by-step case studies showing scider in action. Browse the full list
under the
[**Vignettes**](https://chenlaboratory.github.io/scider/articles/index.html)
tab, or start here:

- [Region based spatial transcriptomics analysis with
  scider](https://chenlaboratory.github.io/scider/articles/getting-started.html):
  density estimation, region detection, and cell-type co-localization on
  a spatial dataset.

## Citation

If you use scider in your work, please cite:

> Li M, Liu N, Nguyen QH, Chen Y (2025). Preserving tissue structure
> through density-based spatial analysis with scider. bioRxiv.
> <https://doi.org/10.1101/2025.09.11.675745>
