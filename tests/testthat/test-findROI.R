test_that("findROI works", {
    data("xenium_bc_spe")
    spe <- gridDensity(spe)
    coi <- c("Breast cancer", "Fibroblasts")
    spe1 <- findROI(spe, coi = coi, method = "walktrap")

    expect_true(all(c("coi", "ngrid.min", "roi") %in% names(spe1@metadata)))
    expect_true(is(spe1@metadata$roi, "DFrame"))
    expect_equal(spe1@metadata$coi, coi)

    spe2 <- findROI(spe,
        coi = coi, method = "walktrap",
        probs = 0.9,
        ngrid.min = 30
    )

    expect_lt(nrow(spe2@metadata$roi), nrow(spe1@metadata$roi))

    expect_error(findROI(spe, coi = coi, method = "xyz"))
    expect_silent(findROI(spe, coi = coi, method = "connected"))

    spe3 <- findROI(spe, coi = coi, diag.nodes = TRUE)

    expect_gt(nrow(spe3@metadata$roi), nrow(spe1@metadata$roi))
})
