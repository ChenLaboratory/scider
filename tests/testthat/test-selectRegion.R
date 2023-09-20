test_that("selectRegion works", {
  data("xenium_bc_spe")
  spe_b <- spe[, SummarizedExperiment::colData(spe)$cell_type == "B cells"]
  dat <- as.data.frame(SpatialExperiment::spatialCoords(spe_b))
  
  expect_silent(selectRegion(dat, x.col = "x_centroid", 
                             y.col = "y_centroid"))

})
