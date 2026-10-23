test_that("positions are converted and overlapping boxes are separated", {
    d <- mini_pruned()
    base <- layout_go_dag(d, verbose = FALSE)
    set.seed(1)
    pos <- data.frame(
        id = base$nodes$id,
        x = base$nodes$x * 96 + stats::runif(nrow(base$nodes), -30, 30),
        y = base$nodes$y * 96 + stats::runif(nrow(base$nodes), -30, 30))
    lay <- layout_from_positions(d, pos, relax = TRUE, verbose = FALSE)
    expect_s3_class(lay, "go_dag_layout")
    expect_equal(nrow(lay$nodes), nrow(d$nodes))
    expect_layout_clean(lay)
    expect_gt(lay$bbox$width, 0)
    expect_gt(lay$bbox$height, 0)
    expect_true(all(vapply(lay$edges$route, nrow, integer(1)) >= 2))
})

test_that("missing positions fall back to the layered layout", {
    d <- mini_pruned()
    pos <- data.frame(id = d$nodes$id[1], x = 10, y = 10)
    lay <- layout_from_positions(d, pos, relax = FALSE, verbose = FALSE)
    expect_equal(nrow(lay$nodes), nrow(d$nodes))
    expect_true(all(is.finite(lay$nodes$x)))
    expect_true(all(is.finite(lay$nodes$y)))
})

test_that("positions can be supplied as a named list or matrix", {
    d <- mini_pruned()
    base <- layout_go_dag(d, verbose = FALSE)
    lst <- stats::setNames(lapply(seq_len(nrow(base$nodes)), function(i) {
        c(x = base$nodes$x[i] * 96, y = base$nodes$y[i] * 96)
    }), base$nodes$id)
    lay1 <- layout_from_positions(d, lst, relax = FALSE, verbose = FALSE)
    mat <- cbind(base$nodes$x * 96, base$nodes$y * 96)
    rownames(mat) <- base$nodes$id
    lay2 <- layout_from_positions(d, mat, relax = FALSE, verbose = FALSE)
    expect_equal(lay1$bbox$width, lay2$bbox$width, tolerance = 1e-8)
})
