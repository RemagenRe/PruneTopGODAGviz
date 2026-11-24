test_that("expand_depth=0 is unchanged and one hop adds direct children", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- mini_dag("BP")
    base <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                         expand_depth = 0, p_cutoff = 1, max_nodes = Inf,
                         verbose = FALSE)
    one <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                        expand_depth = 1, p_cutoff = 1, max_nodes = Inf,
                        verbose = FALSE)
    expect_setequal(base$nodes$id, one$nodes$id[one$nodes$role != "expanded"])
    expect_true("GO:0019684" %in% one$nodes$id)
    expect_false("GO:0009765" %in% one$nodes$id)
    expect_equal(one$nodes$role[one$nodes$id == "GO:0019684"], "expanded")
    f <- tempfile(fileext = ".csv")
    write_prune_table(one, f)
    tab <- utils::read.csv(f)
    expect_equal(tab$action[tab$go_id == "GO:0019684"], "kept_expanded")
})

test_that("expand_depth counts hops per selected term, not global levels", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- mini_dag("BP")
    # GO:0015979 is shallow, GO:0009414 is deeper; one hop each.
    one <- prune_go_dag(e, dag = dag, top_n = 0,
                        keep = c("GO:0015979", "GO:0009414"),
                        expand_depth = 1, p_cutoff = 1, max_nodes = Inf,
                        verbose = FALSE)
    expect_true("GO:0019684" %in% one$nodes$id)   # one hop below 0015979
    expect_false("GO:0009765" %in% one$nodes$id)  # two hops: must not appear
    # A global max-level threshold would have included GO:0009765 here.
    two <- prune_go_dag(e, dag = dag, top_n = 0,
                        keep = c("GO:0015979", "GO:0009414"),
                        expand_depth = 2, p_cutoff = 1, max_nodes = Inf,
                        verbose = FALSE)
    expect_true("GO:0009765" %in% two$nodes$id)
})

test_that("explicit filters win over expand_depth", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- mini_dag("BP")
    dropped <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                            expand_depth = 2, p_cutoff = 1, max_nodes = Inf,
                            drop = "GO:0009765", verbose = FALSE)
    expect_false("GO:0009765" %in% dropped$nodes$id)
    expect_true("GO:0019684" %in% dropped$nodes$id)

    cut <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                        expand_depth = 2, p_cutoff = 1, max_nodes = Inf,
                        drop_subtree = "GO:0019684", verbose = FALSE)
    expect_false("GO:0019684" %in% cut$nodes$id)
    expect_false("GO:0009765" %in% cut$nodes$id)
})

test_that("a removed intermediate term does not leave an orphan descendant", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- mini_dag("BP")
    msgs <- character()
    d <- withCallingHandlers(
        prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                     expand_depth = 2, p_cutoff = 1, max_nodes = Inf,
                     drop = "GO:0019684", verbose = FALSE),
        warning = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_length(msgs, 1L)
    expect_match(msgs, "no longer reach a selected term")
    expect_false("GO:0019684" %in% d$nodes$id)
    # GO:0009765 is two hops below the selected term and only reachable
    # through GO:0019684, so it must not appear as a floating node.
    expect_false("GO:0009765" %in% d$nodes$id)
    expect_equal(d$dropped$action[d$dropped$go_id == "GO:0009765"],
                 "disconnected")
    g <- igraph::graph_from_data_frame(
        d$edges[, c("parent", "child")], directed = TRUE,
        vertices = data.frame(name = d$nodes$id, stringsAsFactors = FALSE))
    expect_equal(igraph::components(g, mode = "weak")$no, 1L)
    # `drop` stays a single-term filter: drop_subtree removes more.
    cut <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                        expand_depth = 2, p_cutoff = 1, max_nodes = Inf,
                        drop_subtree = "GO:0019684", verbose = FALSE)
    expect_false("GO:0009765" %in% cut$nodes$id)
    expect_equal(cut$dropped$action[cut$dropped$go_id == "GO:0009765"],
                 "drop_subtree")
})

test_that("a descendant that keeps another path survives a single-term drop", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- list(
        nodes = data.frame(
            id = c("GO:0015979", "GO:0019684", "GO:0006091", "GO:0009765"),
            name = c("photosynthesis", "light reaction", "precursor",
                     "light harvesting"),
            stringsAsFactors = FALSE),
        edges = data.frame(
            parent = c("GO:0015979", "GO:0019684", "GO:0015979", "GO:0006091"),
            child = c("GO:0019684", "GO:0009765", "GO:0006091", "GO:0009765"),
            type = "is_a", stringsAsFactors = FALSE),
        roots = "GO:0015979", source = "synthetic",
        version = NA_character_, file = NA_character_)
    d <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                      expand_depth = 2, p_cutoff = 1, max_nodes = Inf,
                      drop = "GO:0019684", verbose = FALSE)
    expect_false("GO:0019684" %in% d$nodes$id)
    expect_true("GO:0006091" %in% d$nodes$id)
    expect_true("GO:0009765" %in% d$nodes$id)   # still reachable via GO:0006091
    g <- igraph::graph_from_data_frame(
        d$edges[, c("parent", "child")], directed = TRUE,
        vertices = data.frame(name = d$nodes$id, stringsAsFactors = FALSE))
    expect_equal(igraph::components(g, mode = "weak")$no, 1L)
})

test_that("no single-term drop leaves a disconnected fragment", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- mini_dag("BP")
    full <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                         expand_depth = 3, p_cutoff = 1, max_nodes = Inf,
                         verbose = FALSE)
    for (v in setdiff(full$nodes$id, "GO:0015979")) {
        p <- suppressWarnings(prune_go_dag(
            e, dag = dag, top_n = 0, keep = "GO:0015979", expand_depth = 3,
            p_cutoff = 1, max_nodes = Inf, drop = v, verbose = FALSE))
        expect_true(all(p$nodes$id %in% full$nodes$id), info = v)
        g <- igraph::graph_from_data_frame(
            if (nrow(p$edges)) p$edges[, c("parent", "child")] else
                data.frame(parent = character(), child = character()),
            directed = TRUE,
            vertices = data.frame(name = p$nodes$id, stringsAsFactors = FALSE))
        expect_equal(igraph::components(g, mode = "weak")$no, 1L, info = v)
    }
})

test_that("max_nodes may trim expanded terms but never selected terms", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- mini_dag("BP")
    expect_warning(
        d <- prune_go_dag(e, dag = dag, top_n = 0, keep = "GO:0015979",
                          expand_depth = 2, p_cutoff = 1, max_nodes = 2,
                          verbose = FALSE),
        "could not be reached")
    expect_true("GO:0015979" %in% d$nodes$id)
    expect_true(isTRUE(d$params$max_nodes_exceeded))
})
