test_that("run_go_workflow writes figures and tables", {
    skip_if_not_installed("ragg")
    skip_if_not_installed("svglite")
    td <- tempfile("workflow")
    dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)
    res <- run_go_workflow(
        mini_genes(), mini_annotation(), td, ont = "BP", prefix = "mini",
        obo = system.file("extdata", "mini.obo", package = "PruneTopGODAGviz"),
        verbose = FALSE)
    expect_named(res, "BP")
    expect_true(all(c("pdf", "png") %in%
                        tools::file_ext(res$BP$files)))
    expect_true(file.exists(file.path(td, "mini_BP_enrichment.csv")))
    expect_true(file.exists(file.path(td, "mini_BP_prune_table.csv")))
    expect_true(file.exists(file.path(td, "mini_sessionInfo.txt")))
})
