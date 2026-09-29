test_that("plotCorHeatmap works", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
    spe <- gridDensity(spe, coi = coi)
    spe <- findROI(spe, coi = coi, method = "walktrap")
    set.seed(1)
    result <- corDensity(spe, roi = coi, nperm = 49)
    result_mt <- corDensity(spe, roi = coi, method = "modified.ttest")

    expect_silent(plotCorHeatmap(result$ROI))
    expect_silent(plotCorHeatmap(result$overall))
    # the list returned by corDensity can be passed directly
    expect_silent(plotCorHeatmap(result))
    expect_silent(plotCorHeatmap(result, level = "ROI"))
    expect_silent(plotCorHeatmap(result, roi = "1"))
    expect_silent(plotCorHeatmap(result_mt, stats = "t", level = "ROI"))
    expect_error(plotCorHeatmap(result$overall, roi = "1"))
    expect_error(plotCorHeatmap(list(a = 1)))
    # whole-slide results (roi = NULL) have the same structure
    set.seed(1)
    result0 <- corDensity(spe, coi = coi, nperm = 19)
    expect_silent(plotCorHeatmap(result0))
    expect_silent(plotCorHeatmap(result0, level = "ROI", stats = "shift"))
    expect_silent(plotCorHeatmap(result$ROI, stats = "cor.expected"))
    expect_silent(plotCorHeatmap(result_mt$ROI, stats = "t"))
    expect_error(plotCorHeatmap(result$ROI, stats = "t"))
    expect_silent(plotCorHeatmap(result$overall, stats = "shift"))
    expect_silent(plotCorHeatmap(result$ROI, stats = "shift"))
    expect_error(plotCorHeatmap(result_mt$ROI, stats = "shift"))
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
