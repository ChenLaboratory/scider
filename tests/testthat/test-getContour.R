test_that("getContour works", {
  data("xenium_bc_spe")
  spe <- gridDensity(spe)
  coi <- "Breast cancer"
  spe1 <- getContour(spe, coi = coi)
  
  expect_false("breast_cancer_contour" %in% names(spe@metadata))
  expect_true("breast_cancer_contour" %in% names(spe1@metadata))
  
  expect_error(getContour(spe, coi = "xyz"))
  expect_silent(getContour(spe, coi = coi, bins = 50))
})
