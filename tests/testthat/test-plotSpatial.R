test_that("plotSpatial works", {
    data("xenium_bc_spe")


    expect_silent(plotSpatial(spe,
        group.by = "cell_type",
        pt.shape = ".",
        pt.size = 0.3, pt.alpha = 0.2
    ))
    expect_silent(plotSpatial(spe))

    expect_silent(plotSpatial(spe, reverseY = TRUE))
})
