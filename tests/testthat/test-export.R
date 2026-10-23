test_that("pdf, png, svg and tiff exports are written", {
    skip_if_not_installed("ragg")
    skip_if_not_installed("svglite")
    d <- mini_pruned()
    prefix <- tempfile("go_export")
    files <- save_go_dag(d, prefix, dpi = 100,
                         formats = c("pdf", "png", "svg", "tiff"), bg = "white")
    expect_equal(length(files), 4)
    expect_true(all(file.exists(files)))
    expect_true(all(file.size(files) > 1000))
    expect_true(all(grepl("\\.(pdf|png|svg|tiff)$", files)))
    unlink(files)
})

test_that("a dragged layout can be exported", {
    skip_if_not_installed("ragg")
    d <- mini_pruned()
    base <- layout_go_dag(d, verbose = FALSE)
    pos <- data.frame(id = base$nodes$id, x = base$nodes$x * 96,
                      y = base$nodes$y * 96)
    lay <- layout_from_positions(d, pos, verbose = FALSE)
    p <- plot_go_dag(d, layout = lay)
    expect_s3_class(p, "ggplot")
    f <- file.path(tempdir(), "go_drag")
    files <- save_go_dag(d, f, plot = p, dpi = 100, formats = "png")
    expect_true(file.exists(files[1]))
    unlink(files)
})
