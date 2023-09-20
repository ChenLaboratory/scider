test_that("plotSpatial works", {
    data("xenium_bc_spe")


    expect_silent(plotSpatial(spe,
        shape = ".",
        color = cell_type,
        size = 0.3, alpha = 0.2
    ))
    expect_silent(plotSpatial(spe))

    expect_silent(plotSpatial(spe, reverseY = TRUE))
})
