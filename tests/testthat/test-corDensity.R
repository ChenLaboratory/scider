test_that("corDensity works", {
    data("xenium_bc_spe")

    spe <- gridDensity(spe)

    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")

    spe <- findROI(spe, coi = coi)

    set.seed(1)
    result <- corDensity(spe, roi = coi, nperm = 49)

    expect_true(is(result, "list"))

    expect_true(identical(length(result), 2L))

    expect_true(all(c("celltype1", "celltype2", "cor.coef") %in% names(result$ROI)))

    expect_true(all(c("celltype1", "celltype2", "cor.coef") %in% names(result$overall)))
})

test_that("corDensity permutation method", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
    spe <- gridDensity(spe, coi = coi)
    spe <- findROI(spe, coi = coi, method = "walktrap")

    set.seed(1)
    res1 <- corDensity(spe, roi = coi, nperm = 49)
    set.seed(1)
    res2 <- corDensity(spe, roi = coi, nperm = 49)
    expect_identical(as.data.frame(res1$overall), as.data.frame(res2$overall))

    # the number of threads does not change the results
    set.seed(1)
    res3 <- corDensity(spe, roi = coi, nperm = 49, cpu_threads = 3)
    expect_identical(as.data.frame(res1$overall), as.data.frame(res3$overall))
    expect_identical(as.data.frame(res1$ROI), as.data.frame(res3$ROI))

    cols <-c("celltype1", "celltype2", "cor.coef", "cor.expected", "p.Pos", "p.Neg")
    expect_true(all(cols %in% names(res1$overall)))
    expect_true(all(c("ROI", "ngrid", cols) %in% names(res1$ROI)))
    expect_equal(nrow(res1$overall), 6)

    p <- c(res1$overall$p.Pos, res1$overall$p.Neg)
    expect_true(all(p >= 1 / 50 & p <= 1))
    r <- c(res1$overall$cor.coef, res1$overall$cor.expected)
    expect_true(all(r >= -1 & r <= 1))

    # shift and per pair ROI consistency summary
    ov <- res1$overall
    expect_true(all(c("shift", "n.roi", "n.above", "n.below",
                      "n.sig.Pos", "n.sig.Neg") %in% names(ov)))
    expect_equal(ov$shift, ov$cor.coef - ov$cor.expected)
    nroi <- length(unique(res1$ROI$ROI))
    expect_true(all(ov$n.roi <= nroi))
    expect_true(all(ov$n.above + ov$n.below <= ov$n.roi))
    expect_true(all(ov$n.sig.Pos + ov$n.sig.Neg <= ov$n.roi))
    expect_equal(res1$ROI$shift, res1$ROI$cor.coef - res1$ROI$cor.expected)
    ro <- as.data.frame(res1$ROI)
    ro <- ro[ro$celltype1 == ov$celltype1[1] & ro$celltype2 == ov$celltype2[1], ]
    ok <- !is.na(ro$cor.coef)
    expect_equal(ov$n.above[1], sum(ro$cor.coef[ok] > ro$cor.expected[ok]))
    expect_equal(ov$n.sig.Neg[1], sum(ro$p.Neg[ok] < 0.05))

    # whole-slide path without ROIs
    set.seed(1)
    res0 <- corDensity(spe, coi = coi, nperm = 19)
    expect_named(res0, c("ROI", "overall"))
    expect_identical(unique(res0$ROI$ROI), "1")
    expect_true(all(c("ngrid", "cor.expected", "shift") %in% names(res0$ROI)))
    # a single ROI: the overall results equal the ROI level results
    expect_equal(res0$overall$cor.coef, res0$ROI$cor.coef)
    expect_equal(res0$overall$p.Pos, res0$ROI$p.Pos)
    expect_equal(res0$overall$cor.expected, res0$ROI$cor.expected)

    expect_error(corDensity(spe, roi = coi, nperm = 0), "nperm")

    # grids selected by min.density instead of probs
    ov <- spe@metadata$grid_density$density_overall
    md <- quantile(ov, 0.5)
    set.seed(1)
    r_md <- expect_silent(corDensity(spe, coi = coi, nperm = 9, min.density = md))
    expect_identical(unique(r_md$ROI$ngrid), sum(ov >= md))
    expect_message(corDensity(spe, coi = coi, nperm = 9, probs = 0.5,
                              min.density = md), "min.density")
    expect_error(corDensity(spe, coi = coi, nperm = 9,
                            min.density = max(ov) + 1), "min.density")
    expect_warning(corDensity(spe, roi = coi, nperm = 9, min.density = md),
                   "ignored")
    expect_error(corDensity(spe, roi = coi, id = "not_a_column"), "not found")
})

test_that("corDensity permutation within groups", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
    spe <- gridDensity(spe, coi = coi)
    spe <- findROI(spe, coi = coi, method = "walktrap")

    # a single group is the same as no grouping
    spe$one <- "all"
    set.seed(1); r0 <- corDensity(spe, roi = coi, nperm = 29)
    set.seed(1); r1 <- corDensity(spe, roi = coi, nperm = 29, permute.within = "one")
    expect_identical(as.data.frame(r0$overall), as.data.frame(r1$overall))

    # two groups, each holding two cell types
    spe$compartment <- ifelse(spe$cell_type %in% c("Breast cancer", "Fibroblasts"),
                              "A", "B")
    set.seed(1)
    r2 <- expect_silent(corDensity(spe, roi = coi, nperm = 29, permute.within = "compartment"))
    expect_named(r2, c("ROI", "overall"))
    expect_equal(r2$overall$cor.coef, r0$overall$cor.coef)
    expect_false(isTRUE(all.equal(r2$overall$cor.expected, r0$overall$cor.expected)))

    # groups equal to the cell types: nothing can be permuted
    expect_warning(corDensity(spe, roi = coi, nperm = 9, permute.within = "cell_type"),
                   "never permuted")
    expect_error(corDensity(spe, roi = coi, nperm = 9, permute.within = "nope"), "permute.within")
    spe$bad <- spe$compartment; spe$bad[1] <- NA
    expect_error(corDensity(spe, roi = coi, nperm = 9, permute.within = "bad"), "NA")
    expect_warning(corDensity(spe, roi = coi, method = "modified.ttest",
                              permute.within = "compartment"), "ignored")
})

test_that("corDensity modified.ttest method", {
    data("xenium_bc_spe")
    coi <- c("Breast cancer", "Fibroblasts", "B cells", "T cells")
    spe <- gridDensity(spe, coi = coi)
    spe <- findROI(spe, coi = coi, method = "walktrap")

    res <- corDensity(spe, roi = coi, method = "modified.ttest")
    expect_named(res, c("ROI", "overall"))
    expect_true(all(c("t", "df", "p.Pos", "p.Neg") %in% names(res$ROI)))

    res0 <- corDensity(spe, coi = coi, method = "modified.ttest")
    expect_named(res0, c("ROI", "overall"))
    expect_equal(res0$overall$cor.coef, res0$ROI$cor.coef)
    expect_equal(res0$overall$p.Pos, res0$ROI$p.Pos)
    # one sided p-values recombine to the two sided modified t-test p-value
    expect_equal(pmin(1, 2 * pmin(res$ROI$p.Pos, res$ROI$p.Neg)),
                 2 * pt(-abs(res$ROI$t), res$ROI$df))
})
