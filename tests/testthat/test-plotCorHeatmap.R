test_that("plotCorHeatmap works", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
    spe <- gridDensity(spe, coi = coi)
    spe <- findROI(spe, coi = coi, method = "walktrap")
    result <- corDensity(spe, roi = coi)
    #result <- result$ROI

    expect_silent(plotCorHeatmap(result$ROI))
    expect_silent(plotCorHeatmap(result$overall))
    expect_silent(plotCorHeatmap(result$ROI, stats = "t"))
    expect_error(plotCorHeatmap(result$overall, stats = "xyz"))

    expect_error(plotCorHeatmap(result$ROI, roi = "x"))
    expect_silent(plotCorHeatmap(result$ROI, roi = "1"))
    expect_silent(plotCorHeatmap(result$ROI, roi = c("1", "2", "3")))
    expect_error(plotCorHeatmap(result$ROI, cell.type = "x"))
    expect_silent(plotCorHeatmap(result$ROI,
        cell.type = c(
            "Breast cancer",
            "Fibroblasts"
        )
    ))
})
