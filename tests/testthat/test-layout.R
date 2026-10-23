test_that("the TB layout has no overlaps and only downward edges", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 3, verbose = FALSE)
    lay <- layout_go_dag(d, verbose = FALSE)

    expect_layout_clean(lay)
    lev <- stats::setNames(lay$nodes$level, lay$nodes$id)
    expect_true(all(lev[lay$edges$child] > lev[lay$edges$parent]))
    expect_gt(lay$bbox$width, 0)
    expect_gt(lay$bbox$height, 0)
    expect_true(all(vapply(lay$edges$route, nrow, integer(1)) >= 2))
})

test_that("the LR layout has only rightward edges", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 3, verbose = FALSE)
    lay <- layout_go_dag(d, rankdir = "LR", verbose = FALSE)

    expect_layout_clean(lay)
    x <- stats::setNames(lay$nodes$x, lay$nodes$id)
    expect_true(all(x[lay$edges$child] > x[lay$edges$parent]))
})

test_that("long edges are routed around boxes", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 5, verbose = FALSE)
    lay <- layout_go_dag(d, verbose = FALSE)
    paths <- .edge_paths(lay)
    n <- lay$nodes
    inside <- function(px, py, i) {
        px > n$x[i] - n$w[i] / 2 + 1e-6 && px < n$x[i] + n$w[i] / 2 - 1e-6 &&
            py > -(n$y[i] + n$h[i] / 2) + 1e-6 &&
            py < -(n$y[i] - n$h[i] / 2) - 1e-6
    }
    hits <- 0L
    for (ed in unique(paths$edge)) {
        p <- paths[paths$edge == ed, ]
        ends <- lay$edges[ed, c("parent", "child")]
        for (k in seq_len(nrow(n))) {
            if (n$id[k] %in% ends) next
            if (any(mapply(inside, p$x, p$y, MoreArgs = list(i = k)))) {
                hits <- hits + 1L
                break
            }
        }
    }
    expect_equal(hits, 0)
})
