test_that("spe2PB works", {
    data("xenium_bc_spe")
    spe <- gridDensity(spe)
    coi <- "Breast cancer"
    spe <- findROI(spe, coi = coi)
    spe <- allocateCells(spe)

    y <- spe2PB(spe)

    expect_lt(ncol(y$counts), ncol(spe))

    expect_error(spe2PB(spe, group.id = "xyz"))
})
