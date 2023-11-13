test_that("plotCorHeatmap works", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
    spe <- gridDensity(spe, coi = coi)
    spe <- findROI(spe, coi = coi, method = "walktrap")
    model_result <- corDensity(spe)
    model_result <- model_result$ROI

    expect_silent(plotCorHeatmap(model_result))
    expect_silent(plotCorHeatmap(model_result, stats = "t"))
    expect_error(plotCorHeatmap(model_result, stats = "xyz"))

    expect_error(plotCorHeatmap(model_result, roi = "x"))
    expect_silent(plotCorHeatmap(model_result, roi = "1"))
    expect_silent(plotCorHeatmap(model_result, roi = c("1", "2", "3")))
    expect_error(plotCorHeatmap(model_result, cell.type = "x"))
    expect_silent(plotCorHeatmap(model_result,
        cell.type = c(
            "Breast cancer",
            "Fibroblasts"
        )
    ))
})
