test_that("plotContour works", {
  data("xenium_bc_spe")
  spe <- gridDensity(spe)
  coi <- "Breast cancer"
  spe <- getContour(spe, coi = coi)
  
  expect_message(plotContour(spe, coi = coi))
  
  expect_message(plotContour(spe, coi = coi, size = 0.3, alpha = 0.2))
  
  expect_error(plotContour(spe, coi = "xyz"))
  
  expect_error(plotContour(spe, coi = coi, id = "xyz"))
  
  expect_error(plotContour(spe, coi = coi, overlay = "xyz"))
  
  expect_silent(plotContour(spe, coi = coi, overlay = "density"))
  
  expect_silent(plotContour(spe, coi = coi, overlay = "density",
                            sub.level = "1"))
  
  expect_message(plotContour(spe, coi = coi, sub.level = "1"))
  
  expect_error(plotContour(spe, coi = coi, sub.level = "xyz"))
})
