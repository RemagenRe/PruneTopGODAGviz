test_that("scale_go_layout scales geometry uniformly", {
    d <- mini_pruned()
    lay <- layout_go_dag(d, fontsize = 8, verbose = FALSE)
    lay2 <- scale_go_layout(lay, 2)
    expect_equal(lay2$bbox$width, lay$bbox$width * 2)
    expect_equal(lay2$bbox$height, lay$bbox$height * 2)
    expect_equal(lay2$nodes$x, lay$nodes$x * 2)
    expect_equal(lay2$nodes$w, lay$nodes$w * 2)
    expect_equal(lay2$params$level_gap, lay$params$level_gap * 2)
    expect_error(scale_go_layout(lay, 0), "positive")
})

test_that("fit_go_layout_to_page respects the requested page", {
    d <- mini_pruned()
    fit <- fit_go_layout_to_page(
        d, width_in = 6, height_in = 4, fontsize = 8,
        level_gap = 0.28, sibling_gap = 0.12, edge_style = "bezier")
    expect_s3_class(fit$layout, "go_dag_layout")
    expect_true(is.finite(fit$text_pt) && fit$text_pt > 0)
    expect_lte(fit$size_in["width"], 6 * 1.02)
    expect_lte(fit$size_in["height"], 4 * 1.02)
    expect_true(fit$fits)
})