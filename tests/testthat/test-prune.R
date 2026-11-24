test_that("pruning keeps selected terms and their ancestors only", {
    a <- mini_annotation()
    dag <- mini_dag("BP")
    ego <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(ego, dag = dag, top_n = 3, p_cutoff = 0.05,
                      cluster_jaccard = NULL, verbose = FALSE)

    expect_s3_class(d, "pruned_go_dag")
    expect_true(all(d$nodes$id %in% dag$nodes$id))
    expect_equal(sum(d$nodes$selected), 3)
    expect_true(all(d$nodes$id[d$nodes$selected] %in% ego$result$ID))

    anc <- dag_ancestors(dag, d$nodes$id[d$nodes$selected])
    expect_true(all(anc %in% d$nodes$id))
    expect_true(all(d$edges$parent %in% d$nodes$id))
    expect_true(all(d$edges$child %in% d$nodes$id))
    expect_setequal(d$dropped$go_id, setdiff(dag$nodes$id, d$nodes$id))
    expect_false(any(d$nodes$id == "GO:0000003"))
})

test_that("pruning never changes the statistics", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 4, verbose = FALSE)
    sel <- d$nodes[d$nodes$selected, , drop = FALSE]
    ref <- e$result[match(sel$id, e$result$ID), , drop = FALSE]
    expect_equal(sel$p.adjust, ref$p.adjust)
    expect_equal(sel$pvalue, ref$pvalue)
})

test_that("drop_subtree removes the whole branch", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 5,
                      drop_subtree = "GO:0006950", verbose = FALSE)
    expect_false(any(c("GO:0006950", "GO:0009409", "GO:0009651") %in%
                         d$nodes$id))
    expect_true(all(d$dropped$action[d$dropped$go_id == "GO:0009409"] ==
                        "drop_subtree"))
})

test_that("node cap and chain collapsing work", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 2, max_nodes = 12,
                      collapse_chains = TRUE, cluster_jaccard = NULL,
                      verbose = FALSE)
    expect_lte(nrow(d$nodes), 12)
    expect_true(nrow(d$edges) >= 1)
    expect_true(all(d$edges$n_skipped >= 0))
})

test_that("deleting a parent keeps shared children", {
    dag <- list(
        nodes = data.frame(id = c("root", "A", "B", "C", "D")),
        edges = data.frame(
            parent = c("root", "root", "root", "A", "B", "C"),
            child = c("A", "B", "C", "D", "D", "D"),
            type = "is_a", stringsAsFactors = FALSE))
    expect_setequal(.cascade_orphans(dag, "B"), "B")
    expect_setequal(.cascade_orphans(dag, c("A", "B", "C")),
                    c("A", "B", "C", "D"))
})

test_that("deleting a sole parent removes the child", {
    dag <- list(
        nodes = data.frame(id = c("root", "B", "D")),
        edges = data.frame(parent = c("root", "B"), child = c("B", "D"),
                           type = "is_a", stringsAsFactors = FALSE))
    expect_setequal(.cascade_orphans(dag, "B"), c("B", "D"))
})

test_that("the decision table is written", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    d <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 2, verbose = FALSE)
    f <- tempfile(fileext = ".csv")
    write_prune_table(d, f)
    tab <- utils::read.csv(f)
    expect_true(all(c("go_id", "name", "p.adjust", "action", "reason") %in%
                        names(tab)))
    expect_true(any(tab$action == "kept_selected"))
})
test_that("max_nodes reports node and selected counts separately", {
    anno <- mini_annotation()
    ego <- enrich_go(mini_genes(), anno, ont = "BP", verbose = FALSE)
    msgs <- character()
    d <- withCallingHandlers(
        prune_go_dag(ego, dag = mini_dag("BP"), top_n = 5, p_cutoff = 1,
                     max_nodes = 2, verbose = FALSE),
        warning = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_length(msgs, 1L)
    expect_match(msgs, "could not be reached")
    expect_true(grepl(paste0(nrow(d$nodes), " node(s) remain"), msgs,
                      fixed = TRUE))
    expect_true(grepl(paste0(sum(d$nodes$selected),
                             " selected term(s) protected"), msgs,
                      fixed = TRUE))
    expect_true(sum(d$nodes$selected) < nrow(d$nodes))
    expect_gt(nrow(d$nodes), 2)
    expect_true(isTRUE(d$params$max_nodes_exceeded))
})

test_that("direct and propagated counts match their definitions", {
    a <- mini_annotation()
    e <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    dag <- mini_dag("BP")
    args <- list(x = e, dag = dag, top_n = 0, keep = "GO:0015979",
                 expand_depth = 2, p_cutoff = 1, max_nodes = Inf,
                 verbose = FALSE)
    direct <- do.call(prune_go_dag, c(args, list(counts = "direct")))
    prop <- do.call(prune_go_dag, c(args, list(counts = "propagated")))
    sets <- split(a$pairs$transcript, a$pairs$go_id)
    for (i in seq_len(nrow(direct$nodes))) {
        id <- direct$nodes$id[[i]]
        own <- if (is.null(sets[[id]])) character() else sets[[id]]
        expect_equal(direct$nodes$n_in_term[[i]],
                     length(intersect(own, e$gene)), info = id)
        expect_equal(direct$nodes$n_bg_term[[i]],
                     length(intersect(own, e$universe)), info = id)
        ids <- dag_descendants(dag, id)   # true-path rule, term included
        expect_true(id %in% ids, info = id)
        pool <- unique(unlist(sets[ids], use.names = FALSE))
        expect_equal(prop$nodes$n_in_term[[i]],
                     length(intersect(pool, e$gene)), info = id)
        expect_equal(prop$nodes$n_bg_term[[i]],
                     length(intersect(pool, e$universe)), info = id)
    }
    expect_true(any(prop$nodes$n_in_term > direct$nodes$n_in_term))
    expect_equal(direct$nodes$n_in_universe[[1]], length(e$gene))
    expect_equal(direct$nodes$n_bg_universe[[1]], length(e$universe))
})