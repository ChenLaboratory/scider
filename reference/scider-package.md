# scider: Spatial cell-type inter-correlation by density in R

scider is an user-friendly R package providing functions to model the
global density of cells in a slide of spatial transcriptomics data. All
functions in the package are built based on the SpatialExperiment
object, allowing integration into various spatial
transcriptomics-related packages from Bioconductor. After modelling
density, the package allows for several downstream analysis, including
colocalization analysis, boundary detection analysis and differential
density analysis.

`scider` implements functions to analyse spatial transcriptomics data
with cell type annotations by performing cell type correlation via
density estimation and cell type co-localization via real number
distance. Functions include density estimation, statistical modelling
and visualizations.

## Details

`scider` uses SpatialExperiment objects as the main infrastructure,
which can easily be integrated with a wide variety of Bioconductor
packages.

## See also

Useful links:

- <https://github.com/ChenLaboratory/scider>

- <https://chenlaboratory.github.io/scider/>

- Report bugs at <https://github.com/ChenLaboratory/scider/issues>

## Author

**Maintainer**: Yunshun Chen <yuchen@wehi.edu.au>
([ORCID](https://orcid.org/0000-0003-4911-5653))

Authors:

- Mengbo Li <li.me@wehi.edu.au>
  ([ORCID](https://orcid.org/0000-0002-9666-5810))

- Ning Liu <liu.n@wehi.edu.au>
  ([ORCID](https://orcid.org/0000-0002-9487-9305))

- Quoc Hoang Nguyen <nguyen.q@wehi.edu.au>
  ([ORCID](https://orcid.org/0009-0007-2828-7567))

Ning Liu <liu.n@wehi.edu.au>, Mengbo Li <li.me@wehi.edu.au>, Yunshun
Chen <yuchen@wehi.edu.au>, Quoc Hoang Nguyen <nguyen.q@wehi.edu.au>
