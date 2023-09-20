test_that("gridDensity works", {
  data("xenium_bc_spe")
  
  spe1 <- gridDensity(spe)
  
  expect_equal(length(spe@metadata),0)
  
  expect_gt(length(spe1@metadata),0)
  
  expect_true("grid_density" %in% names(spe1@metadata))
  
  expect_true("grid_info" %in% names(spe1@metadata))
  
  expect_silent(gridDensity(spe, coi = c("Breast cancer", "Fibroblasts")))
  
  expect_error(gridDensity(spe, coi = c("breast cancer", "fibroblast")))
  
  expect_error(gridDensity(spe, coi = c("Breast cancer", "Fibroblasts"), 
                           id = "cell_typess"))
  
  expect_error(gridDensity(spe, kernel = "xyz"))
  
  expect_silent(gridDensity(spe, bandwidth = 100))
})
