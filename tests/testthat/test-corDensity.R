test_that("corDensity works", {
    data("xenium_bc_spe")

    spe <- gridDensity(spe)

    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")

    spe <- findROI(spe, coi = coi)

    result <- corDensity(spe)

    expect_true(is(result, "DFrame"))

    expect_true(all(c("celltype1", "celltype2", "cor.coef") %in% names(result)))

    result2 <- corDensity(spe, by.roi = FALSE)

    expect_true(all(c("celltype1", "celltype2", "cor.coef") %in% names(result2)))

    expect_lt(nrow(result2), nrow(result))
})
