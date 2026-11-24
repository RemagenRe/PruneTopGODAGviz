test_that("bezier export ignores intermediate layout waypoints", {
    route <- rbind(c(0, 0), c(0, 1), c(1, 2))
    lay <- list(
        edges = data.frame(parent = "A", child = "D", type = "is_a",
                           n_skipped = 0L, stringsAsFactors = FALSE),
        edge_style = "bezier", rankdir = "TB")
    lay$edges$route <- list(route)
    lay$edges$route_ortho <- list(route)
    lay2 <- lay
    lay2$edges$route <- list(route[c(1L, 3L), , drop = FALSE])

    p1 <- PruneTopGODAGviz:::.edge_paths(lay)
    p2 <- PruneTopGODAGviz:::.edge_paths(lay2)
    expect_equal(p1[, c("x", "y")], p2[, c("x", "y")])
    expect_equal(unique(p1$edge), 1L)
})

test_that("layout_from_positions keeps the requested edge style", {
    d <- mini_pruned()
    pos <- data.frame(id = d$nodes$id,
                      x = seq_len(nrow(d$nodes)) * 100,
                      y = seq_len(nrow(d$nodes)) * 50)
    lay <- layout_from_positions(d, pos, edge_style = "bezier",
                                 relax = FALSE, verbose = FALSE)
    expect_equal(lay$edge_style, "bezier")
})
