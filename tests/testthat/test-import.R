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