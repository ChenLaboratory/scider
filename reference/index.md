# Package index

## All functions

- [`allocateCells()`](https://chenlaboratory.github.io/scider/reference/allocateCells.md)
  : Annotate all cells with contour level of cell type-specific density.
- [`cellsInRegion()`](https://chenlaboratory.github.io/scider/reference/cellsInRegion.md)
  : Check which cells are in which regions
- [`computeDensity()`](https://chenlaboratory.github.io/scider/reference/computeDensity.md)
  : Perform kernel density estimation on SpatialExperiment
- [`computeDensityHex()`](https://chenlaboratory.github.io/scider/reference/computeDensityHex.md)
  : Perform kernel density estimation on SpatialExperiment
- [`contour2sf()`](https://chenlaboratory.github.io/scider/reference/contour2sf.md)
  : Draw a contour region on some density level
- [`coord_hash()`](https://chenlaboratory.github.io/scider/reference/coord_hash.md)
  : Hash two 15-bytes signed integers into one 32-bytes integer.
- [`corDensity()`](https://chenlaboratory.github.io/scider/reference/corDensity.md)
  : Test for density correlation between two cell types.
- [`findNbrsGrid()`](https://chenlaboratory.github.io/scider/reference/findNbrsGrid.md)
  : Construct a neighbour list from grid coordinates.
- [`findNbrsSNN()`](https://chenlaboratory.github.io/scider/reference/findNbrsSNN.md)
  : Construct a SNN neighbour list from assay.
- [`findNbrsSpatial()`](https://chenlaboratory.github.io/scider/reference/findNbrsSpatial.md)
  : Construct a distance-based neighbour list from cell coordinates.
- [`findROI()`](https://chenlaboratory.github.io/scider/reference/findROI.md)
  : Find ROIs based on cell type-specific densities via graph-based
  method.
- [`getClusters()`](https://chenlaboratory.github.io/scider/reference/getClusters.md)
  : Cluster cells in spe using graph methods.
- [`getContour()`](https://chenlaboratory.github.io/scider/reference/getContour.md)
  : Get contour from density
- [`getContourRegions()`](https://chenlaboratory.github.io/scider/reference/getContourRegions.md)
  : Calculate areas between every two density levels
- [`getDE()`](https://chenlaboratory.github.io/scider/reference/getDE.md)
  : Differential expression between two clusters or two groups of
  clusters.
- [`getHVG()`](https://chenlaboratory.github.io/scider/reference/getHVG.md)
  : Get top highly variable genes.
- [`getMarkers()`](https://chenlaboratory.github.io/scider/reference/getMarkers.md)
  : Find up-regulated marker genes for all clusters (quasi-NB GLM).
- [`getNiche()`](https://chenlaboratory.github.io/scider/reference/getNiche.md)
  : Build a niche assay based on the profile of neighbouring cells
- [`getSubClusters()`](https://chenlaboratory.github.io/scider/reference/getSubClusters.md)
  : Sub-cluster a single cluster.
- [`globalMoran()`](https://chenlaboratory.github.io/scider/reference/globalMoran.md)
  : Calculate global Moran for 1 to 2 variables.
- [`grid2df()`](https://chenlaboratory.github.io/scider/reference/grid2df.md)
  : Convert x,y nodes to data.frame of polygons
- [`grid2sf()`](https://chenlaboratory.github.io/scider/reference/grid2sf.md)
  : Convert x,y nodes to sf polygons
- [`gridDensity()`](https://chenlaboratory.github.io/scider/reference/gridDensity.md)
  : Perform kernel density estimation on SpatialExperiment for cell
  types of interest
- [`gridSPE()`](https://chenlaboratory.github.io/scider/reference/gridSPE.md)
  : Summarize a SpatialExperiment object at grid-level
- [`localMoran()`](https://chenlaboratory.github.io/scider/reference/localMoran.md)
  : Calculate local Moran for 1 to 2 variables.
- [`mergeClusters()`](https://chenlaboratory.github.io/scider/reference/mergeClusters.md)
  : Manually merge clusters.
- [`mergeROI()`](https://chenlaboratory.github.io/scider/reference/mergeROI.md)
  : Manually merge ROIs
- [`normalizeAssay()`](https://chenlaboratory.github.io/scider/reference/normalizeAssay.md)
  : Perform log normalization for counts
- [`plotCellCompo()`](https://chenlaboratory.github.io/scider/reference/plotCellCompo.md)
  : Plot cell type composition in each density level of cell of
  interest.
- [`plotContour()`](https://chenlaboratory.github.io/scider/reference/plotContour.md)
  : Plot contour lines.
- [`plotContourRegion()`](https://chenlaboratory.github.io/scider/reference/plotContourRegion.md)
  : Visualising an sf object (for internal use only at the moment)
- [`plotCorHeatmap()`](https://chenlaboratory.github.io/scider/reference/plotCorHeatmap.md)
  : Plot model statistics using heatmap.
- [`plotDR()`](https://chenlaboratory.github.io/scider/reference/plotDR.md)
  [`plotUMAP()`](https://chenlaboratory.github.io/scider/reference/plotDR.md)
  [`plotPCA()`](https://chenlaboratory.github.io/scider/reference/plotDR.md)
  : Plot reduced dimensions.
- [`plotDensCor()`](https://chenlaboratory.github.io/scider/reference/plotDensCor.md)
  : Plot density correlation between two cell types
- [`plotDensity()`](https://chenlaboratory.github.io/scider/reference/plotDensity.md)
  : Plot grid-based density.
- [`plotDots()`](https://chenlaboratory.github.io/scider/reference/plotDots.md)
  : Dot plot of gene expression by groups
- [`plotGrid()`](https://chenlaboratory.github.io/scider/reference/plotGrid.md)
  : Plot grid from metadata.
- [`plotImage()`](https://chenlaboratory.github.io/scider/reference/plotImage.md)
  : Plot background image of spe
- [`plotLISA()`](https://chenlaboratory.github.io/scider/reference/plotLISA.md)
  : Plotting LISA (e.g. moran)
- [`plotLISAscatter()`](https://chenlaboratory.github.io/scider/reference/plotLISAscatter.md)
  : Scatterplot for local moran's I
- [`plotROI()`](https://chenlaboratory.github.io/scider/reference/plotROI.md)
  : Plot ROIs on spatial.
- [`plotSpatial()`](https://chenlaboratory.github.io/scider/reference/plotSpatial.md)
  : Plot cells based on spatial coordinates.
- [`plotTopMarkers()`](https://chenlaboratory.github.io/scider/reference/plotTopMarkers.md)
  : Dot plot or heatmap of the top marker genes per cluster.
- [`plotViolin()`](https://chenlaboratory.github.io/scider/reference/plotViolin.md)
  : Violin plot using genes or cell data
- [`postSelRegion()`](https://chenlaboratory.github.io/scider/reference/postSelRegion.md)
  : Merge sel_region from the selectRegion function to
  SpatialExperiment.
- [`readProseg()`](https://chenlaboratory.github.io/scider/reference/readProseg.md)
  : Read Proseg V2 output into spe
- [`readVisium()`](https://chenlaboratory.github.io/scider/reference/readVisium.md)
  : Read Visium output into spe
- [`readVisiumHD()`](https://chenlaboratory.github.io/scider/reference/readVisiumHD.md)
  : Read VisiumHD output into spe
- [`readXenium()`](https://chenlaboratory.github.io/scider/reference/readXenium.md)
  : Read Xenium output into spe
- [`realignVisium()`](https://chenlaboratory.github.io/scider/reference/realignVisium.md)
  : Scale and straighten out Visium coordinates
- [`realignVisiumHD()`](https://chenlaboratory.github.io/scider/reference/realignVisiumHD.md)
  : Scale and straighten out VisiumHD coordinates
- [`relabelClusters()`](https://chenlaboratory.github.io/scider/reference/relabelClusters.md)
  : Renumber all clusters by size.
- [`runPCA()`](https://chenlaboratory.github.io/scider/reference/runPCA.md)
  : Fast PCA using irlba.
- [`runUMAP()`](https://chenlaboratory.github.io/scider/reference/runUMAP.md)
  : UMAP using uwot. Parameters are set to be similar to Seurat's
- [`selectRegion()`](https://chenlaboratory.github.io/scider/reference/selectRegion.md)
  : Select region of interest from plot
- [`spe2PB()`](https://chenlaboratory.github.io/scider/reference/spe2PB.md)
  : Given a 'SpatialExperiment' data object, create pseudo-bulk samples
  using the colData information and return a DGEList object
- [`` `[`( ``*`<SpatialExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://chenlaboratory.github.io/scider/reference/sub-SpatialExperiment-ANY-ANY-ANY-method.md)
  : Subset for grid level analysis
- [`topMarkers()`](https://chenlaboratory.github.io/scider/reference/topMarkers.md)
  : Combine the top markers of each cluster into one data frame.
- [`update_bound()`](https://chenlaboratory.github.io/scider/reference/update_bound.md)
  : Update the x,y limits of a plot
