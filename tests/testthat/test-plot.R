test_that("plot_go_dag returns a ggplot with a computed size", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 3, verbose = FALSE)
    p <- plot_go_dag(d)
    expect_s3_class(p, "ggplot")
    sz <- attr(p, "go_size")
    expect_true(is.list(sz))
    expect_gt(sz$width, 1)
    expect_gt(sz$height, 1)
})

test_that("save_go_dag writes pdf and png files", {
    skip_if_not_installed("ragg")
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 2, verbose = FALSE)
    prefix <- file.path(tempdir(), "go_test")
    files <- save_go_dag(d, prefix, dpi = 300)
    expect_equal(length(files), 2)
    expect_true(all(file.exists(files)))
    expect_true(all(file.size(files) > 1000))
    unlink(files)
})

test_that("interactive output can be created", {
    skip_if_not_installed("ggiraph")
    skip_if_not_installed("htmlwidgets")
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 2, verbose = FALSE)
    f <- file.path(tempdir(), "go_test.html")
    w <- if (nzchar(Sys.which("pandoc"))) {
        plot_go_dag_interactive(d, file = f)
    } else {
        suppressWarnings(plot_go_dag_interactive(d, file = f))
    }
    expect_s3_class(w, "girafe")
    expect_true(file.exists(f))
    unlink(f)
})
