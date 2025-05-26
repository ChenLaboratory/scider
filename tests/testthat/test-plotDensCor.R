test_that("plotDensCor works", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts")
    spe <- gridDensity(spe, coi = coi)
    spe <- findROI(spe, coi = coi, method = "walktrap")

    expect_silent(plotDensCor(spe,
        celltype1 = "Breast cancer",
        celltype2 = "Fibroblasts"
    ))
    expect_error(plotDensCor(spe,
        celltype1 = "xyz",
        celltype2 = "Fibroblasts"
    ))
    expect_silent(plotDensCor(spe,
        celltype1 = "Breast cancer",
        celltype2 = "Fibroblasts", roi = coi
    ))
    expect_silent(plotDensCor(spe,
        celltype1 = "Breast cancer",
        celltype2 = "Fibroblasts", fit = "linear"
    ))
    expect_silent(plotDensCor(spe,
        celltype1 = "Breast cancer",
        celltype2 = "Fibroblasts", df = 2
    ))
})
