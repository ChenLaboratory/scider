test_that("plotROI works", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts")
    spe <- gridDensity(spe, coi = coi)
    expect_error(plotROI(spe))

    spe <- findROI(spe, coi = coi, method = "walktrap", steps = 5)

    expect_silent(plotROI(spe, roi = coi, pt.size = 0.3, pt.alpha = 0.2))
    expect_silent(plotROI(spe, roi = coi))
    expect_silent(plotROI(spe, roi = coi, show.legend = TRUE))
    #expect_error(plotROI(spe, roi = coi, id = "xyz"))
})
