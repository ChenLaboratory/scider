# scider 1.8.0 (2025-10-30)
* Update help pages
* New function findNbrsGrid() to construct a neighbour list from grid coordinates
* New function findNbrsSNN() to construct a SNN neighbour list from assay
* New function findNbrsSpatial() to construct a distance-based neighbour list from cell coordinates
* New function getNiche() to build a niche assay based on the profile of neighbouring cells
* New function getClusters() to cluster cells in spe using graph methods
* New functions globalMoran() and localMoran() to calculate Moran's I statistics
* New function normalizeAssay() to perform spatial data normalization
* New function plotDR() for plotting reduced dimensions
* New functions plotImage() and plotLISA() for visualization
* New function realignVisium() to scale and straighten out Visium coordinates
* New functions runPCA() and runUMAP() for performing dimensionality reduction

# scider 1.4.0 (2024-10-30)
* gridDensity() now supports hexagonal grids.
* Make 'grid.length.x' = 100 by default in gridDensity() and remove the default value for 'ngrid.x'.
* New functions plotGrid() and gridSPE()
* gridSPE() has an option to calculate grid-level count matrix for each cell type
* Remove redundant argument 'n' in plotSpatial()
* Rename 'coi' to 'contour', 'all' to 'self.included' in plotCellCompo()
* Fixed a bug in plotContour() where 'overlay' is not working properly
* Fixed a bug in corDensity() that produced NaN values
* New argument 'contour' in allocateCells() allowing user to specify the contour levels for cell allocation. Default to all detected contours.
* corDensity() can be run without running findROI(), and it returns a DataFrame instead of a list of DataFrames if no ROI detected.
* New argument 'probs' added to plotDensity() for grid filtering when ROIs are not identified previously.
* 'equal.cell' is set to TRUE by default in getContour().
* plotContour() now plots contour for all the cells when 'coi' is NULL (by default).
* When 'coi' contains multiple cell types, plotContour() looks for the corresponding contour in spe@metadata for plotting.
* Added another option 'none' for 'overlay' in plotContour().
* For multiple 'coi', getContour() now includes all coi in the contour labelling.
* gridDensity() returns an extra column 'density_overall' in spe@metadata$grid_density.
* plotDensity() plots overall cell density by default.
* The default value of 'probs' in plotDensity() is now 0.5.
* For multiple 'coi', plotDensity() now plots the total instead of average.
* Fixed a ggtitle() bug in plotDensity() when 'coi' contains more than one cell type.
* findROI() now use all cells (overall) by default.
* coord_fixed for both plotDensity() and plotROI().
* plotROI() now plots all the cells in the background when 'roi' is 'overall'.
* In corDensity(), 'trace' is now set to FALSE by default.
* Fixed a bug in cellsInRegion() when there is only one ROI.
* Fixed bugs in contour2sf() when no 'areas_up'.
* Converted cell numbers as a vector in spe2PB().
* Fixed a bug in spe2PB() when the count table from spe is not matrix.
* Fixed a bug in plotContour() to ensure exact match for cell types to avoid similar cell type names.

# scider 1.2.0 (2023-10-30)
* Bioconductor release.

# scider 0.99.0
* Initial submission to _Bioconductor_.