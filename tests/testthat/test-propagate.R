test_that("propagation applies the true-path rule", {
    a <- mini_annotation()
    d <- mini_dag("BP")
    p <- propagate_go_annotation(a, d, ont = "BP", verbose = FALSE)
    expect_s3_class(p, "go_annotation")
    expect_true(isTRUE(p$propagated))
    bp <- a$pairs[a$pairs$go_id %in% a$terms$go_id[a$terms$ontology == "BP"], ]
    expect_gt(nrow(p$pairs), nrow(bp))
    expect_true(all(c("GO:0008150", "GO:0008152") %in% p$terms$go_id))
    expect_true("GO:0008150" %in% p$pairs$go_id)
    expect_equal(p$universe$BP, a$universe$BP)
})

test_that("propagation keeps only terms present in the DAG", {
    a <- mini_annotation()
    d <- mini_dag("BP")
    p <- propagate_go_annotation(a, d, ont = "BP", verbose = FALSE)
    expect_true(all(p$pairs$go_id %in% d$nodes$id))
    expect_true(all(p$terms$go_id %in% d$nodes$id))
})
