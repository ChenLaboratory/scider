test_that("corDensity works", {
    data("xenium_bc_spe")

    spe <- gridDensity(spe)

    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")

    spe <- findROI(spe, coi = coi)

    result <- corDensity(spe)
    
    expect_true(is(result, "list"))
    
    expect_true(identical(length(result), 2L))

    expect_true(all(c("celltype1", "celltype2", "cor.coef") %in% names(result$ROI)))
    
    expect_true(all(c("celltype1", "celltype2", "cor.coef") %in% names(result$overall)))
})
