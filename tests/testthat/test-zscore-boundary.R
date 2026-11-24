test_that("zScore is NA when the hypergeometric variance is zero", {
    t2g <- data.frame(go_id = "T1", transcript = "g1", stringsAsFactors = FALSE)
    t2n <- data.frame(go_id = "T1", name = "term one", stringsAsFactors = FALSE)
    res <- PruneTopGODAGviz:::.ora_enrich("g1", "g1", t2g, t2n, "BH", 1, 10)
    expect_true(is.na(res$result$zScore[[1]]))
    expect_true(is.finite(res$result$pvalue[[1]]))

    t2g2 <- data.frame(go_id = "T1", transcript = c("g1", "g2"),
                       stringsAsFactors = FALSE)
    res2 <- PruneTopGODAGviz:::.ora_enrich(c("g1", "g2"), c("g1", "g2"),
                                            t2g2, t2n, "BH", 1, 10)
    expect_true(is.na(res2$result$zScore[[1]]))
    expect_true(is.finite(res2$result$pvalue[[1]]))
})

test_that("zScore stays finite for a normal background", {
    t2g <- data.frame(
        go_id = c(rep("T1", 10), rep("T2", 10)),
        transcript = c(sprintf("g%02d", 1:10), sprintf("g%02d", 11:20)),
        stringsAsFactors = FALSE)
    t2n <- data.frame(go_id = c("T1", "T2"), name = c("term one", "term two"),
                      stringsAsFactors = FALSE)
    res <- PruneTopGODAGviz:::.ora_enrich(sprintf("g%02d", 1:5),
                                          sprintf("g%02d", 1:20),
                                          t2g, t2n, "BH", 1, 30)
    expect_true(is.finite(res$result$zScore[[1]]))
})
