test_that("expanded non-selected leaves can be merged by jaccard", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- mini_dag("BP")
    d <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0006950",
                      expand_depth = 1, cluster_jaccard = 0.8,
                      p_cutoff = 1, max_nodes = Inf, verbose = FALSE)
    kids <- c("GO:0009409", "GO:0009651")
    expect_equal(sum(kids %in% d$nodes$id), 1L)
    expect_true("GO:0006950" %in% d$nodes$id)
    expect_true(any(d$nodes$n_collapsed >= 1L))
})

test_that("jaccard merging is skipped when no annotation is available", {
    dag <- mini_dag("BP")
    ids <- dag$nodes$id[seq_len(5)]
    tab <- data.frame(
        Term = ids, Description = paste("term", seq_along(ids)),
        PValue = seq(0.01, 0.05, length.out = length(ids)),
        FDR = seq(0.02, 0.06, length.out = length(ids)),
        Count = seq_along(ids), BgCount = 10 + seq_along(ids),
        stringsAsFactors = FALSE)
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    y <- import_go_enrich(f, ontology = "BP", dag = dag, verbose = FALSE)
    expect_warning(
        d <- prune_go_dag(y, dag = dag, top_n = 5, p_cutoff = 1,
                          cluster_jaccard = 0.8, verbose = FALSE),
        "no annotation")
    expect_s3_class(d, "pruned_go_dag")
})
