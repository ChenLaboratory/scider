test_that("multiplication works", {
    data("xenium_bc_spe")
    spe <- gridDensity(spe)
    spe <- findROI(spe, coi = c("Breast cancer", "Fibroblasts"))
    spe <- getContour(spe, coi = "Breast cancer")
    spe <- allocateCells(spe)

    expect_silent(plotCellCompo(spe, coi = "Breast cancer"))

    expect_silent(plotCellCompo(spe, coi = "Breast cancer", by.roi = TRUE))

    expect_error(plotCellCompo(spe, coi = "xyz"))

    expect_error(plotCellCompo(spe, coi = "Breast cancer", id = "xyz"))

    expect_silent(plotCellCompo(spe,
        coi = "Breast cancer",
        level.name = "breast_cancer_contour"
    ))

    expect_error(plotCellCompo(spe,
        coi = "Breast cancer",
        level.name = "xyz_contour"
    ))
})
