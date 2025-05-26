test_that("allocateCells works", {
    data("xenium_bc_spe")
    spe <- gridDensity(spe)
    coi <- "Breast cancer"
    spe <- findROI(spe, coi = coi)

    spe1 <- allocateCells(spe)

    expect_gt(ncol(spe1@colData), ncol(spe@colData))

    expect_true("breast_cancer_roi" %in% names(spe1@colData))

    spe <- getContour(spe, coi)

    spe2 <- allocateCells(spe)

    expect_gt(ncol(spe2@colData), ncol(spe1@colData))

    expect_false("breast_cancer_contour" %in% names(spe1@colData))

    expect_true("breast_cancer_contour" %in% names(spe2@colData))
})
