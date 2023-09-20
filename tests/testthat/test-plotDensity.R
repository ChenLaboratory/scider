test_that("plotDensity works", {
  data("xenium_bc_spe")
  spe <- gridDensity(spe)
  
  expect_silent(plotDensity(spe, coi = "Breast cancer"))
  expect_silent(plotDensity(spe, coi = "Fibroblasts"))
  expect_error(plotDensity(spe, coi = "xxx"))
  expect_silent(plotDensity(spe, coi = "Fibroblasts", probs = 0.5))
  expect_error(plotDensity(spe, coi = "Fibroblasts", probs = "xx"))
})
