test_that("jaccard merging collapses a highly similar sibling leaf", {
    nodes <- data.frame(id = c("A", "L1", "L2", "L3"),
                        name = c("parent", "one", "two", "three"),
                        selected = FALSE, stringsAsFactors = FALSE)
    edges <- data.frame(parent = c("A", "A", "A"),
                        child = c("L1", "L2", "L3"),
                        type = "is_a", stringsAsFactors = FALSE)
    annotation <- list(pairs = data.frame(
        go_id = c("L1", "L1", "L2", "L2", "L3"),
        transcript = c("g1", "g2", "g1", "g2", "g3"),
        stringsAsFactors = FALSE))
    stat <- c(A = 1, L1 = 0.2, L2 = 0.1, L3 = 0.3)
    z <- .cluster_jaccard(nodes, edges, selected = character(),
                          annotation = annotation, stat_map = stat,
                          threshold = 0.9)
    expect_equal(length(z$merged), 1L)
    expect_true("L1" %in% z$nodes$id)
    expect_false("L2" %in% z$nodes$id)
    expect_true(sum(z$n_collapsed) >= 1L)
})

test_that("layout_from_positions fills missing and duplicate coordinates", {
    d <- mini_pruned()
    ids <- d$nodes$id
    pos <- data.frame(id = ids[1:2], x = c(10, NA), y = c(20, 30))
    lay <- layout_from_positions(d, pos, relax = FALSE, verbose = FALSE)
    expect_equal(nrow(lay$nodes), nrow(d$nodes))
    expect_true(all(is.finite(lay$nodes$x)))
    expect_true(all(is.finite(lay$nodes$y)))

    dup <- data.frame(id = rep(ids[[1]], 2), x = c(1, 2), y = c(3, 4))
    lay2 <- layout_from_positions(d, dup, relax = FALSE, verbose = FALSE)
    expect_equal(nrow(lay2$nodes), nrow(d$nodes))
    expect_true(all(is.finite(lay2$nodes$x)))
    expect_true(all(is.finite(lay2$nodes$y)))
})
