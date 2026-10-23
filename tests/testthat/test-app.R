test_that("the app UI is built in English with Arial", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("visNetwork")
    skip_if_not_installed("colourpicker")
    td <- tempfile("builtin-data")
    dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)
    old <- getOption("PruneTopGODAGviz.data")
    on.exit(options(PruneTopGODAGviz.data = old), add = TRUE)
    file.create(file.path(td, c(
        "At_GO_info_panther.csv", "Maize_GO_info_panther.csv",
        "Rice_GO_info_panther.csv", "Soybean_GO_info_panther.csv",
        "Homo_sapiens_GO_info_panther.csv",
        "Mus_musculus_GO_info_panther.csv")))
    options(PruneTopGODAGviz.data = td)
    ui <- app_ui()
    expect_s3_class(ui, "shiny.tag.list")
    html <- if (requireNamespace("htmltools", quietly = TRUE)) {
        rt <- htmltools::renderTags(ui); paste(rt$head, rt$html)
    } else {
        paste(as.character(ui), collapse = "\n")
    }
    expect_true(grepl("PruneTopGODAGviz Studio", html))
    expect_true(grepl("Run Enrichment", html))
    expect_true(grepl("Delete Term and Descendants", html))
    expect_true(grepl("Enrichment Results (CSV)", html, fixed = TRUE))
    expect_true(grepl("Pruning Decisions (CSV)", html, fixed = TRUE))
    expect_true(grepl("Line Spacing", html, fixed = TRUE))
    expect_true(grepl("Or Use a Built-in Annotation", html, fixed = TRUE))
    expect_true(grepl("Download Built-in GO Info", html, fixed = TRUE))
    expect_true(grepl("arabidopsis_thaliana", html, fixed = TRUE))
    expect_true(grepl("zea_mays", html, fixed = TRUE))
    expect_true(grepl("oryza_sativa", html, fixed = TRUE))
    expect_true(grepl("glycine_max", html, fixed = TRUE))
    expect_true(grepl("homo_sapiens", html, fixed = TRUE))
    expect_true(grepl("mus_musculus", html, fixed = TRUE))
    expect_true(grepl("download_human", html, fixed = TRUE))
    expect_true(grepl("download_mouse", html, fixed = TRUE))
    expect_true(grepl("Propagate Annotations to Ancestors", html, fixed = TRUE))
    expect_true(grepl("#E9C46A", html, fixed = TRUE))
    expect_true(grepl("#EA8379", html, fixed = TRUE))
    expect_false(grepl("At_GO_info.csv", html, fixed = TRUE))
    expect_false(grepl("GO_info_example.png", html, fixed = TRUE))
    expect_true(grepl("Arial", html))
})

test_that("the interactive graph carries events, labels and colours", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("visNetwork")
    d <- mini_pruned()
    v <- .build_vis(
        d, physics = TRUE, spring = 150, gravity = -60, fontsize = 7,
        colors = list(low = "#CFEAF1", mid = "#F27970", high = "#c82423",
                      ancestor = "#FFFFFF", border = "#000000",
                      edge = "#000000", part_of = "#2878B5", bg = "#FFFFFF"))
    expect_s3_class(v, "visNetwork")
    expect_equal(nrow(v$x$nodes), nrow(d$nodes))
    expect_equal(nrow(v$x$edges), nrow(d$edges))
    expect_true(all(c("click", "doubleClick", "dragEnd") %in%
                        names(v$x$events)))
    expect_true(any(grepl("\\n", v$x$nodes$label)))
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", v$x$nodes$color.background)))
    expect_true(any(v$x$edges$dashes))
    expect_true(all(c("x", "y") %in% names(v$x$nodes)))
    expect_true(all(is.finite(v$x$nodes$x)) && all(is.finite(v$x$nodes$y)))
})

test_that("labels follow the show switches", {
    d <- mini_pruned()
    lab <- .vis_labels(d$nodes, 7, show_name = FALSE, show_counts = FALSE,
                       show_p = FALSE)
    expect_true(all(grepl("^GO:[0-9]{7}\\n$", lab)))
})

test_that("the colour scale defaults to FDR 0.05 / 0.01 / 0.0001", {
    expect_equal(.fdr_position(0.05), 0)
    expect_equal(.fdr_position(0.01), 0.5)
    expect_equal(.fdr_position(0.0001), 1)
    expect_equal(.fdr_position(0.00001), 1)
    expect_equal(.fdr_position(0.5), 0)
    nodes <- data.frame(selected = rep(TRUE, 5),
                        p.adjust = c(0.0001, 0.00001, 0.01, 0.05, 0.50))
    cols <- .vis_colours(
        nodes, list(low = "#2878B5", mid = "#FFFFFF",
                    high = "#c82423", ancestor = "#CFEAF1"),
        fdr_breaks = c(0.05, 0.01, 0.0001))
    expect_identical(cols[1], cols[2])
    expect_identical(toupper(cols[3]), "#FFFFFF")
    expect_identical(toupper(cols[4]), "#2878B5")
    expect_identical(toupper(cols[5]), "#2878B5")
})

test_that("the upload limit allows annotation files up to 2 GB", {
    old <- getOption("shiny.maxRequestSize")
    on.exit(options(shiny.maxRequestSize = old), add = TRUE)
    .set_upload_limit()
    expect_equal(getOption("shiny.maxRequestSize"), 2 * 1024^3 - 1)
})

test_that("enrichment results can be exported as CSV", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    f <- tempfile(fileext = ".csv")
    utils::write.csv(e$result, f, row.names = FALSE)
    tab <- utils::read.csv(f)
    expect_true(all(c("ID", "Description", "GeneRatio", "pvalue",
                      "p.adjust", "Count") %in% names(tab)))
    expect_gt(nrow(tab), 0)
})

test_that("statistics are formatted without padding spaces", {
    expect_identical(.fmt_p(3.6e-8), "3.6e-08")
    lab <- paste0("GO:0000001 (", .fmt_p(3.6e-8), ")")
    expect_identical(lab, "GO:0000001 (3.6e-08)")
    expect_false(grepl("\\( ", lab))
})

test_that("run_app checks its suggested packages", {
    expect_true(is.function(run_app))
    expect_true(all(c("port", "launch.browser") %in% names(formals(run_app))))
})

test_that("ID matching deduplicates, sorts and rejects sparse overlaps", {
    ids <- paste0("gene", 1:10)
    x <- .check_id_match(c(ids, ids, "  "), ids)
    expect_equal(x$ids, sort(ids))
    expect_equal(x$n, 10)
    expect_equal(x$n_matched, 10)
    expect_equal(x$rate, 1)

    expect_error(.check_id_match(ids[1:3], ids),
                 class = "id_check_error")
    expect_error(.check_id_match(ids, paste0("other", 1:5)),
                 class = "id_check_error")
    expect_error(.check_id_match(ids, paste0("gene", 1:2)),
                 class = "id_check_error")
})
