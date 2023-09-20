test_that("postSelRegion works", {
  data("xenium_bc_spe")
  coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
  
  spe <- gridDensity(spe, coi = coi)
  
  sel_region <- data.frame(
     "node" = seq(10),
     "node_x" = seq(10),
     "node_y" = seq(10)
  )
  spe1 <- postSelRegion(spe, sel_region)
  
  expect_true("user_defined" %in% spe1@metadata$roi$component)
})
