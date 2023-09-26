# scider: Spatial cell-type inter-correlation by density in R. <img src="man/figures/scider_sticker.png" align="right" alt="" width="120" />

[![R-CMD-check](https://github.com/ChenLaboratory/scider/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/ChenLaboratory/scider/actions)
[![Codecov test coverage](https://codecov.io/gh/ChenLaboratory/scider/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/ChenLaboratory/scider?branch=devel)

*scider* implements functions to analyse spatial transcriptomics 
data with cell type annotations by performing cell type 
correlation via density estimation and cell type co-localization 
via real number distance. Functions include density 
estimation, statistical modelling and visualizations.

Install released version form Bioconductor

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scider")
```


Install development version from GitHub

```
library(devtools)   
devtools::install_github("ChenLaboratory/scider")
```
