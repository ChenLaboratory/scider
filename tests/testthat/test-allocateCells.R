test_that("allocateCells works", {
    data("xenium_bc_spe")
    spe <- gridDensity(spe)
    coi <- "Breast cancer"
    spe <- findROI(spe, coi = coi)

    spe1 <- allocateCells(spe)

    expect_gt(ncol(colData(spe1)), ncol(colData(spe)))

    expect_true("roi" %in% names(colData(spe1)))

    spe <- getContour(spe, coi)

    spe2 <- allocateCells(spe)

    expect_gt(ncol(colData(spe2)), ncol(colData(spe1)))

    expect_false("breast_cancer_contour" %in% names(colData(spe1)))

    expect_true("breast_cancer_contour" %in% names(colData(spe2)))
})
