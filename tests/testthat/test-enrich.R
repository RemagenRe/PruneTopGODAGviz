test_that("enrich_go uses a hypergeometric over-representation test", {
    a <- mini_annotation()
    g <- mini_genes()
    ego <- enrich_go(g, a, ont = "BP", verbose = FALSE)

    terms <- a$terms[a$terms$ontology == "BP", , drop = FALSE]
    pairs <- a$pairs[a$pairs$go_id %in% terms$go_id, , drop = FALSE]
    hit <- intersect(g, a$universe$BP)
    uni <- a$universe$BP
    t2g <- unique(pairs[, c("go_id", "transcript")])
    sets <- split(t2g$transcript, t2g$go_id)
    q <- unique(t2g[t2g$transcript %in% hit, , drop = FALSE])
    q_terms <- unique(q$go_id)
    sets <- sets[q_terms]
    sets <- lapply(sets, intersect, uni)
    keep <- lengths(sets) >= 10 & lengths(sets) <= 500
    sets <- sets[keep]
    q_terms <- names(sets)
    k <- vapply(q_terms, function(x) {
        length(intersect(q$transcript[q$go_id == x], uni))
    }, integer(1))
    M <- lengths(sets)
    N <- length(uni)
    n <- length(hit)
    p <- stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    expect_equal(sort(ego$result$pvalue), unname(sort(p)))
    expect_equal(sort(ego$result$p.adjust), unname(sort(stats::p.adjust(p, "BH"))))
    expect_equal(ego$result$GeneRatio,
                 sprintf("%s/%s", k[match(ego$result$ID, q_terms)], n))
    expect_equal(ego$result$BgRatio,
                 sprintf("%s/%s", M[match(ego$result$ID, q_terms)], N))
    expect_equal(ego$result$Count, unname(k[match(ego$result$ID, q_terms)]))
    expect_true(all(c("FoldEnrichment", "RichFactor", "BgCount") %in%
                        names(ego$result)))
    expect_false(anyNA(ego$result$qvalue))
    expect_true(all(ego$result$qvalue >= 0 & ego$result$qvalue <= 1))
})

test_that("unannotated input identifiers trigger a warning", {
    a <- mini_annotation()
    expect_warning(enrich_go(c(mini_genes(), "NOT_A_GENE"), a, ont = "BP",
                             verbose = FALSE), "not annotated")
    expect_error(enrich_go(rep("NOT_A_GENE", 5), a, ont = "BP",
                           verbose = FALSE), "None of the input")
})

test_that("as_go_enrich wraps a result data frame", {
    a <- mini_annotation()
    ego <- enrich_go(mini_genes(), a, ont = "BP", verbose = FALSE)
    wrapped <- as_go_enrich(ego$result, a, ont = "BP")
    expect_s3_class(wrapped, "go_enrich")
    expect_equal(nrow(wrapped$result), nrow(ego$result))
    expect_true(all(c("FoldEnrichment", "BgCount") %in% names(wrapped$result)))
})
test_that("custom universe excludes unannotated query IDs from GeneRatio", {
    a <- mini_annotation()
    g <- mini_genes()
    uni <- c(a[["universe"]][["BP"]], "NOT_A_GENE")
    e <- enrich_go(c(g, "NOT_A_GENE"), a, ont = "BP", universe = uni,
                   verbose = FALSE)
    den <- as.numeric(sub(".*/", "", e$result$GeneRatio))
    mapped <- length(intersect(g, a[["universe"]][["BP"]]))
    expect_true(all(den == mapped))
    expect_false("NOT_A_GENE" %in% e$gene)
    p <- prune_go_dag(e, dag = mini_dag("BP"), top_n = 5, p_cutoff = 1,
                      verbose = FALSE)
    expect_true(all(p$nodes$n_in_universe == mapped))
    expect_true(all(p$nodes$n_bg_universe == length(e$universe)))
})