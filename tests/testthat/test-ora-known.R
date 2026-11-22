test_that("ORA matches a hand-computed hypergeometric example", {
    universe <- sprintf("g%03d", 1:100)
    pairs <- rbind(
        data.frame(go_id = "T1", transcript = universe[1:20]),
        data.frame(go_id = "T3", transcript = universe[1:5]),
        data.frame(go_id = "T0", transcript = universe[21:100]),
        stringsAsFactors = FALSE)
    terms <- data.frame(
        go_id = c("T1", "T3", "T0"),
        name = c("term one", "term three", "background filler"),
        ontology = "BP", stringsAsFactors = FALSE)
    anno <- list(pairs = pairs, terms = terms,
                 universe = list(BP = universe),
                 keytype = "transcript",
                 source = list(file = "synthetic", size = NA, mtime = NA,
                               n_rows = nrow(pairs), n_bad = 0L, sep = ","))
    class(anno) <- "go_annotation"

    hit <- universe[1:10]
    ego <- enrich_go(hit, anno, ont = "BP", minGSSize = 1,
                     maxGSSize = 1000, verbose = FALSE)
    res <- ego$result

    p1 <- stats::phyper(10 - 1, 20, 80, 10, lower.tail = FALSE)
    p3 <- stats::phyper(5 - 1, 5, 95, 10, lower.tail = FALSE)
    expect_equal(nrow(res), 2L)
    expect_false("T0" %in% res$ID)  # zero hits: not tested
    expect_equal(res$pvalue[res$ID == "T1"], p1)
    expect_equal(res$pvalue[res$ID == "T3"], p3)
    expect_equal(res$p.adjust,
                 unname(stats::p.adjust(c(p1, p3), "BH")[
                     match(res$ID, c("T1", "T3"))]))
    expect_equal(res$GeneRatio[res$ID == "T1"], "10/10")
    expect_equal(res$GeneRatio[res$ID == "T3"], "5/10")
    expect_equal(res$BgRatio[res$ID == "T1"], "20/100")
    expect_equal(res$BgRatio[res$ID == "T3"], "5/100")
    expect_equal(res$Count[res$ID == "T1"], 10L)
    expect_equal(res$Count[res$ID == "T3"], 5L)
    expect_equal(res$BgCount[res$ID == "T1"], 20L)
    expect_equal(res$BgCount[res$ID == "T3"], 5L)
    expect_equal(res$FoldEnrichment[res$ID == "T1"], 5)
    expect_equal(res$FoldEnrichment[res$ID == "T3"], 10)
    expect_equal(ego$params$n_terms_tested, 2L)
    expect_equal(ego$params$test_scope, "hit_terms")
})

test_that("imported FDR is passed through and never recalculated", {
    tab <- data.frame(
        ID = c("GO:0006412", "GO:0006508"),
        Description = c("translation", "proteolysis"),
        pvalue = c(0.001, 0.02),
        FDR = c(0.0123, 0.4),
        Count = c(12, 4), BgCount = c(300, 50),
        stringsAsFactors = FALSE)
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    x <- parse_external_table(f, ontology = "BP")
    expect_equal(x$result$p.adjust, c(0.0123, 0.4))
    expect_equal(x$result$qvalue, x$result$p.adjust)
    expect_equal(x$pvalue_kind, "source_adjusted")
})
