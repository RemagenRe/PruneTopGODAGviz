test_that("PANTHER parser handles mixed group result shapes", {
    entry <- function(id, label, n, nr, p, fdr) {
        list(
            input_list = list(
                number_in_list = n, pValue = p, fdr = fdr,
                fold_enrichment = 2, expected = 1,
                mapped_id_list = list(mapped_id = c("AT1G00010", "AT1G00020"))
            ),
            number_in_reference = nr,
            term = list(level = 0, id = id, label = label)
        )
    }
    obj <- list(overrepresentation = list(
        annotation_type = "GO biological process complete",
        group = list(
            list(result = list(
                entry("GO:0006412", "translation", 2, 300, 1e-8, 1e-6),
                entry("GO:0006508", "proteolysis", 2, 50, 1e-4, 1e-3)
            )),
            list(result = entry("GO:0006396", "RNA processing", 2, 100, 1e-5, 1e-4))
        )
    ))
    f <- tempfile(fileext = ".json")
    jsonlite::write_json(obj, f, auto_unbox = TRUE)
    x <- parse_panther_json(f)
    expect_equal(nrow(x$result), 3)
    expect_equal(x$ontology, "BP")
    expect_true(all(c("ID", "Description", "pvalue", "p.adjust", "qvalue",
                      "Count", "BgCount", "geneID") %in% names(x$result)))
    expect_equal(x$result$qvalue, x$result$p.adjust)
    expect_equal(x$result$count_label[1], "2 | 300")
    expect_equal(x$result$p.adjust[1], 1e-6)
})

test_that("generic external table accepts common column names", {
    tab <- data.frame(
        Term = c("GO:0006412", "GO:0006508"),
        Description = c("translation", "proteolysis"),
        PValue = c(1e-8, 1e-4),
        FDR = c(1e-6, 1e-3),
        Count = c(12, 4),
        BgCount = c(300, 50),
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    x <- parse_external_table(f, ontology = "BP")
    expect_equal(nrow(x$result), 2)
    # Mixed-case headers must map to the right columns, not to each other.
    expect_equal(x$result$ID, c("GO:0006412", "GO:0006508"))
    expect_equal(x$result$pvalue, c(1e-8, 1e-4))
    expect_equal(x$result$p.adjust, c(1e-6, 1e-3))
    expect_equal(x$result$Count, c(12, 4))
    expect_equal(x$result$BgCount, c(300, 50))
    expect_equal(x$result$p.adjust, x$result$qvalue)
    expect_equal(x$counts$count_label, c("12 | 300", "4 | 50"))
})
test_that("external counts use raw Count | BgCount in pruning", {
    dag <- mini_dag("BP")
    ids <- dag$nodes$id[seq_len(2)]
    res <- data.frame(
        ID = ids,
        Description = c("term one", "term two"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.02, 0.03),
        qvalue = c(0.02, 0.03),
        Count = c(12, 4),
        BgCount = c(300, 50),
        geneID = c("A/B", "C"),
        stringsAsFactors = FALSE)
    ext <- list(
        result = res, enrich = NULL, gene = c("A", "B", "C"),
        universe = character(), ontology = "BP", annotation = NULL,
        params = list(keytype = "external"),
        external_counts = data.frame(
            ID = ids, n_in_term = res$Count, n_bg_term = res$BgCount,
            n_in_universe = NA_real_, n_bg_universe = NA_real_,
            count_label = c("12 | 300", "4 | 50"),
            stringsAsFactors = FALSE))
    class(ext) <- "go_enrich"
    p <- prune_go_dag(ext, dag = dag, top_n = 2, p_cutoff = 1,
                      verbose = FALSE)
    kept <- p$nodes$id %in% ids
    expect_equal(p$nodes$label_l3[kept],
                 c("12 | 300", "4 | 50"))
    expect_true(all(p$nodes$label_l3[!kept] == ""))
})
test_that("generic table ontology is read from a namespace column", {
    tab <- data.frame(
        Term = c("GO:0004672", "GO:0005524"),
        Description = c("protein kinase activity", "ATP binding"),
        Namespace = c("molecular_function", "molecular_function"),
        PValue = c(1e-9, 1e-6),
        FDR = c(2e-7, 4e-5),
        Count = c(25, 60),
        BgCount = c(900, 1500),
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    x <- parse_external_table(f)
    expect_equal(x$ontology, "MF")
})

test_that("import stops when no term exists in the chosen ontology", {
    tab <- data.frame(
        Term = "GO:0004672", Description = "protein kinase activity",
        PValue = 1e-9, FDR = 2e-7, Count = 25, BgCount = 900,
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    dag <- mini_dag("BP")
    expect_error(
        import_go_enrich(f, ontology = "BP", dag = dag, verbose = FALSE),
        "Choose the ontology")
})

test_that("ontology text helper maps common spellings", {
    expect_equal(PruneTopGODAGviz:::.ontology_from_text("molecular_function"), "MF")
    expect_equal(PruneTopGODAGviz:::.ontology_from_text("cellular component"), "CC")
    expect_equal(PruneTopGODAGviz:::.ontology_from_text("BP"), "BP")
})
test_that("external import enforces the file size limit", {
    expect_equal(PruneTopGODAGviz:::.max_external_bytes, 50 * 1024^2)
    dag <- mini_dag("BP")
    ids <- dag$nodes$id[seq_len(2)]
    tab <- data.frame(
        Term = ids, Description = c("term one", "term two"),
        PValue = c(1e-9, 1e-6), FDR = c(2e-7, 4e-5),
        Count = c(25, 60), BgCount = c(900, 1500),
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    expect_error(
        import_go_enrich(f, ontology = "BP", dag = dag,
                         max_bytes = 10, verbose = FALSE),
        "limit is")
    expect_s3_class(
        import_go_enrich(f, ontology = "BP", dag = dag,
                         max_bytes = 1e6, verbose = FALSE),
        "go_enrich")
})
test_that("ratio columns contribute the numerator, not the denominator", {
    dag <- mini_dag("BP")
    ids <- dag$nodes$id[seq_len(2)]
    tab <- data.frame(
        Term = ids, Description = c("term one", "term two"),
        PValue = c(1e-4, 0.02), FDR = c(1e-3, 0.04),
        GeneRatio = c("12/100", "4/80"),
        BgRatio = c("300/20000", "50/5000"),
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    x <- parse_external_table(f, ontology = "BP")
    expect_equal(x$result$Count, c(12, 4))
    expect_equal(x$result$BgCount, c(300, 50))
    expect_equal(x$result$RichFactor, c(12 / 300, 4 / 50))
    expect_equal(x$counts$count_label, c("12 | 300", "4 | 50"))
})

test_that("Count that disagrees with the ratio numerator is reported", {
    dag <- mini_dag("BP")
    tab <- data.frame(
        Term = dag$nodes$id[[1]], Description = "term one",
        PValue = 1e-4, FDR = 1e-3, Count = 12,
        GeneRatio = "10/100", BgCount = 300, BgRatio = "300/20000",
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    expect_warning(parse_external_table(f, ontology = "BP"),
                   "differs from the GeneRatio numerator")
    x <- suppressWarnings(parse_external_table(f, ontology = "BP"))
    expect_equal(x$result$Count, 12)
})

test_that("malformed ratios become NA with a warning", {
    dag <- mini_dag("BP")
    tab <- data.frame(
        Term = dag$nodes$id[seq_len(3)], Description = c("a", "b", "c"),
        PValue = c(0.01, 0.02, 0.03), FDR = c(0.04, 0.05, 0.06),
        GeneRatio = c("12/0", "a/b", ""),
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    expect_warning(parse_external_table(f, ontology = "BP"),
                   "invalid GeneRatio")
    x <- suppressWarnings(parse_external_table(f, ontology = "BP"))
    expect_true(all(is.na(x$result$Count)))
    expect_true(all(!nzchar(x$result$count_label)))
})

test_that("a raw p-value never impersonates an adjusted p-value", {
    dag <- mini_dag("BP")
    ids <- dag$nodes$id[seq_len(2)]
    tab <- data.frame(
        Term = ids, Description = c("term one", "term two"),
        PValue = c(0.03, 0.5), Count = c(12, 4), BgCount = c(300, 50),
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    x <- parse_external_table(f, ontology = "BP")
    expect_equal(x$result$pvalue, c(0.03, 0.5))
    expect_true(all(is.na(x$result$p.adjust)))
    expect_true(all(is.na(x$result$qvalue)))
    expect_equal(x$pvalue_kind, "raw_only")

    y <- import_go_enrich(f, ontology = "BP", dag = dag, verbose = FALSE)
    expect_equal(y$params$stat_default, "pvalue")
    p <- prune_go_dag(y, dag = dag, stat = "pvalue", top_n = 2,
                      p_cutoff = 0.05, verbose = FALSE)
    expect_true(any(p$nodes$selected))
    expect_true(all(is.na(p$nodes$p.adjust[p$nodes$selected])))
})

test_that("missing FDR entries stay missing", {
    dag <- mini_dag("BP")
    tab <- data.frame(
        Term = dag$nodes$id[seq_len(2)], Description = c("a", "b"),
        PValue = c(0.001, 0.02), FDR = c(0.01, NA),
        Count = c(12, 4), BgCount = c(300, 50),
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    x <- parse_external_table(f, ontology = "BP")
    expect_equal(x$result$p.adjust, c(0.01, NA))
    expect_equal(x$pvalue_kind, "source_adjusted_with_missing")
    expect_equal(x$n_missing_adjusted, 1L)
})
test_that("import and pruning record provenance and the used statistic", {
    dag <- mini_dag("BP")
    ids <- dag$nodes$id[seq_len(2)]
    tab <- data.frame(
        Term = ids, Description = c("term one", "term two"),
        PValue = c(1e-3, 0.02), FDR = c(1e-2, 0.05),
        Count = c(12, 4), BgCount = c(300, 50),
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".csv")
    utils::write.csv(tab, f, row.names = FALSE)
    y <- import_go_enrich(f, ontology = "BP", dag = dag, verbose = FALSE)
    expect_true(nzchar(y$params$package_version))
    expect_true(nzchar(y$params$r_version))
    p <- prune_go_dag(y, dag = dag, stat = "p.adjust", top_n = 2,
                      p_cutoff = 1, verbose = FALSE)
    expect_equal(p$provenance$stat, "p.adjust")
    expect_true(nzchar(p$provenance$package_version))
    expect_equal(p$nodes$color_value, p$nodes$p.adjust)
})