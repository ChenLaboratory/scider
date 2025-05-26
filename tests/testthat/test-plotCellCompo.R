test_that("multiplication works", {
    data("xenium_bc_spe")
    spe <- gridDensity(spe)
    spe <- findROI(spe, coi = c("Breast cancer", "Fibroblasts"))
    spe <- getContour(spe, coi = "Breast cancer")
    spe <- allocateCells(spe)

    expect_silent(plotCellCompo(spe, contour = "Breast cancer"))

    expect_silent(plotCellCompo(spe, contour = "Breast cancer", roi = c("Breast cancer", "Fibroblasts")))

    expect_error(plotCellCompo(spe, contour = "xyz"))

    expect_error(plotCellCompo(spe, contour = "Breast cancer", id = "xyz"))


    # expect_silent(plotCellCompo(spe,
    #     contour = "Breast cancer",
    #     level.name = "breast_cancer_contour"
    # ))
    # 
    # expect_error(plotCellCompo(spe,
    #     contour = "Breast cancer",
    #     level.name = "xyz_contour"
    # ))
})
