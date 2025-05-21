test_that("mergeROI works", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts")
    spe <- gridDensity(spe, coi = coi)
    spe <- findROI(spe, coi = coi, method = "walktrap")
    spe1 <- mergeROI(spe, list("1-2" = 1:2))

    expect_false("component_before_merge" %in% names(spe@metadata$roi))
    expect_true("component_before_merge" %in% names(spe1@metadata$roi))
    expect_false("1-2" %in% spe@metadata$roi$component)
    expect_true("1-2" %in% spe1@metadata$roi$component)
})
